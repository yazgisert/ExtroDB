# Flags the species among the non-extremophile proteins that meet one of our extremophile thresholds, based on the
# growth optima of their type strain in BacDive and on TEMPURA.

import csv
import json
import re
import time
from collections import Counter
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.parse import quote

from bacdive_client import BACDIVE_API_BASE, ensure_list, fetch_json
from curate_halophile_bacdive import parse_concentration_to_molar

NONEXTREM_FAA = "nonextremophile_data/nonextremophile_filtered.faa"
TEMPURA_CSV   = "curation/source_tables/200617_TEMPURA.csv"
CACHE_JSONL   = "curation/source_tables/bacdive_nonextremophile_cache.jsonl"
OUT_TSV       = "curation/source_tables/nonextremophile_species_flags.tsv"

REQUEST_DELAY = 0.1
BATCH_SIZE    = 100

# Thresholds of the curation (Table 3.1); the acid threshold is the highest optimum among our acidotolerant organisms
THERMOPHILE_MIN_OGT   = 45.0
PSYCHROPHILE_MAX_OGT  = 15.0
HALOPHILE_MIN_NACL    = 0.2
ACIDOPHILE_MAX_PH     = 6.0
ALKALIPHILE_MIN_PH    = 9.0
ALKALITOLERANT_MIN_PH = 7.5   # optimum above 7.5 and documented growth at pH >= 9.0


def species_name(organism):
    # "Thermus thermophilus (strain HB8)" -> "Thermus thermophilus"
    return " ".join(re.sub(r"\(.*?\)", "", organism).replace('"', "").split()[:2])


def read_species():
    counts = Counter()
    with open(NONEXTREM_FAA) as f:
        for line in f:
            if line.startswith(">"):
                match = re.search(r" OS=(.+?) OX=", line)
                if match:
                    counts[species_name(match.group(1))] += 1
    return counts


def read_tempura():
    temperatures = {}
    with open(TEMPURA_CSV) as f:
        for row in csv.DictReader(f):
            try:
                temperatures.setdefault(species_name(row["genus_and_species"]), []).append(float(row["Topt_ave"]))
            except ValueError:
                continue
    return {name: sum(values) / len(values) for name, values in temperatures.items()}


def midpoint(text):
    # "28", "5.5-6.5" or "22 - 30" -> mean of the numbers
    numbers = [float(x) for x in re.findall(r"\d+(?:\.\d+)?", str(text))]
    return sum(numbers) / len(numbers) if numbers else None


def get(url):
    for attempt in range(3):
        try:
            time.sleep(REQUEST_DELAY)
            return fetch_json(url)
        except HTTPError as exc:
            if exc.code == 404:
                return None
        except (URLError, TimeoutError, ConnectionError):
            pass
        time.sleep(5 * (attempt + 1))
    raise SystemExit(f"BacDive request failed three times: {url}. Run the script again to resume.")


def bacdive_ids(species):
    genus, epithet = species.split()
    ids, url = [], f"{BACDIVE_API_BASE}/taxon/{quote(genus)}/{quote(epithet)}"
    while url:
        payload = get(url)
        if not payload:
            break
        ids += [str(i) for i in payload.get("results", [])]
        url = payload.get("next")
    return ids


def strain_conditions(record):
    culture = record.get("Culture and growth conditions", {})
    type_strain = str(record.get("Name and taxonomic classification", {}).get("type strain", "")).lower() == "yes"
    conditions = {"temperature_optimum": [], "ph_optimum": [], "ph_growth": [], "nacl_optimum": []}
    for entry in ensure_list(culture.get("culture temp")):
        if isinstance(entry, dict) and "opt" in str(entry.get("type", "")).lower():
            conditions["temperature_optimum"].append(midpoint(entry.get("temperature")))
    for entry in ensure_list(culture.get("culture pH")):
        if not isinstance(entry, dict):
            continue
        if "opt" in str(entry.get("type", "")).lower():
            conditions["ph_optimum"].append(midpoint(entry.get("pH")))
        elif str(entry.get("ability", "")).lower() == "positive":
            numbers = [float(x) for x in re.findall(r"\d+(?:\.\d+)?", str(entry.get("pH")))]
            conditions["ph_growth"] += numbers
    for entry in ensure_list(record.get("Physiology and metabolism", {}).get("halophily")):
        if (isinstance(entry, dict) and str(entry.get("salt", "")).lower() == "nacl"
                and str(entry.get("tested relation", "")).lower() == "optimum"):
            concentration = parse_concentration_to_molar(str(entry.get("concentration", "")))
            if concentration:
                conditions["nacl_optimum"].append(concentration["midpoint_M"])
    return {"type_strain": type_strain, **{k: [v for v in values if v is not None] for k, values in conditions.items()}}


def fetch_strains(species):
    ids, strains = bacdive_ids(species), []
    for start in range(0, len(ids), BATCH_SIZE):
        payload = get(f"{BACDIVE_API_BASE}/fetch/{';'.join(ids[start:start + BATCH_SIZE])}")
        if payload:
            strains += [strain_conditions(record) for record in payload.get("results", {}).values()
                        if isinstance(record, dict)]
    return strains


def flags(strains, tempura_ogt, type_strain_only=True):
    found = set()
    for s in strains:
        if type_strain_only and not s.get("type_strain"):
            continue
        if any(t >= THERMOPHILE_MIN_OGT for t in s["temperature_optimum"]):
            found.add("thermophile")
        if any(t <= PSYCHROPHILE_MAX_OGT for t in s["temperature_optimum"]):
            found.add("psychrophile")
        if any(c >= HALOPHILE_MIN_NACL for c in s["nacl_optimum"]):
            found.add("halophile")
        if any(p <= ACIDOPHILE_MAX_PH for p in s["ph_optimum"]):
            found.add("acidophile")
        if any(p >= ALKALIPHILE_MIN_PH for p in s["ph_optimum"]):
            found.add("alkaliphile")
        elif (any(p > ALKALITOLERANT_MIN_PH for p in s["ph_optimum"])
              and any(p >= ALKALIPHILE_MIN_PH for p in s["ph_growth"])):
            found.add("alkaliphile")
    if tempura_ogt is not None and tempura_ogt >= THERMOPHILE_MIN_OGT:
        found.add("thermophile")
    if tempura_ogt is not None and tempura_ogt <= PSYCHROPHILE_MAX_OGT:
        found.add("psychrophile")
    return sorted(found)


def summary(values):
    return ";".join(f"{v:g}" for v in sorted(set(round(v, 2) for v in values)))


if __name__ == "__main__":
    species = read_species()
    tempura = read_tempura()
    cache = {}
    if Path(CACHE_JSONL).exists():
        with open(CACHE_JSONL) as f:
            for line in f:
                entry = json.loads(line)
                cache[entry["species"]] = entry["strains"]
    # Cache entries from before the type strain was recorded are queried again, but only if some strain meets a
    # threshold; otherwise the type strain cannot meet one either
    todo = [s for s in species if len(s.split()) == 2 and (
        s not in cache or (any("type_strain" not in st for st in cache[s])
                           and flags(cache[s], None, type_strain_only=False)))]
    print(f"{len(species):,} species, {len(cache):,} already in the cache, {len(todo):,} to query")

    with open(CACHE_JSONL, "a") as f:
        for i, name in enumerate(todo, 1):
            cache[name] = fetch_strains(name)
            f.write(json.dumps({"species": name, "strains": cache[name]}) + "\n")
            f.flush()
            if i % 50 == 0:
                print(f"  {i:,} / {len(todo):,}")

    rows = []
    for name, n_proteins in species.most_common():
        strains = cache.get(name, [])
        type_strains = [s for s in strains if s.get("type_strain")]
        merged = {k: [v for s in type_strains for v in s[k]] for k in ("temperature_optimum", "ph_optimum", "nacl_optimum")}
        rows.append({
            "species": name, "n_proteins": n_proteins, "n_bacdive_strains": len(strains),
            "n_type_strains": len(type_strains),
            "temperature_optima": summary(merged["temperature_optimum"]), "ph_optima": summary(merged["ph_optimum"]),
            "nacl_optima_M": summary(merged["nacl_optimum"]),
            "tempura_ogt": f"{tempura[name]:g}" if name in tempura else "",
            "flags": ",".join(flags(strains, tempura.get(name))),
        })
    with open(OUT_TSV, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    in_bacdive = [r for r in rows if r["n_bacdive_strains"]]
    flagged = [r for r in rows if r["flags"]]
    print(f"\n{len(in_bacdive):,} of {len(rows):,} species found in BacDive")
    print(f"{len(flagged):,} species flagged, {sum(r['n_proteins'] for r in flagged):,} of "
          f"{sum(species.values()):,} proteins")
    for cls in ("thermophile", "psychrophile", "halophile", "acidophile", "alkaliphile"):
        hits = [r for r in flagged if cls in r["flags"]]
        print(f"  {cls:<13} {len(hits):>4} species  {sum(r['n_proteins'] for r in hits):>6} proteins")
    print(f"\nSaved to {OUT_TSV}")
