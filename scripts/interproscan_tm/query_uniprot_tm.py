# Fetches transmembrane annotations from UniProt: extremophile WP_ IDs are first mapped to UniProtKB,
# non-extremophile accessions are queried directly. Resumable via a progress file.

import json
import random
import time
import urllib.parse
import urllib.request
import urllib.error
from collections import defaultdict
from pathlib import Path

SEED          = 42
N_NONEXTREM   = 500          # non-extremophile sequences to sample
BATCH_SIZE    = 50           # IDs per UniProt search query
MAP_BATCH     = 500          # IDs per ID-mapping job
POLL_INTERVAL = 10           # seconds between ID-mapping polls
MAX_POLL      = 60
MAX_RETRIES   = 5

random.seed(SEED)

BASE_UNIPROT  = "https://rest.uniprot.org"
OUT_DIR       = Path("results/tm_analysis")
PROGRESS_F    = OUT_DIR / "progress.json"

PROK_SAMPLE   = Path("results/interproscan/prokaryote_sample.faa")
NONEXT_FAA    = Path("nonextremophile_data/nonextremophile_filtered.faa")

CLASS_NAMES   = ["acidophile", "alkaliphile", "halophile", "psychrophile", "thermophile"]


def read_fasta_ids(path):
    records = []
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                parts = line[1:].strip().split("||")
                pid   = parts[0].split()[0]
                label = parts[-1].strip() if len(parts) > 1 else ""
                records.append((pid, label))
    return records


def read_nonextrem_sample(faa_path, n=N_NONEXTREM):
    all_ids = []
    with open(faa_path) as f:
        for line in f:
            if line.startswith(">"):
                pid = line[1:].strip().split("||")[0].split()[0]
                all_ids.append(pid)
    rng = random.Random(SEED)
    sampled = rng.sample(all_ids, min(n, len(all_ids)))
    print(f"  Sampled {len(sampled):,} non-extremophile IDs (from {len(all_ids):,})")
    return sampled


def load_progress():
    if PROGRESS_F.exists():
        with open(PROGRESS_F) as f:
            return json.load(f)
    return {
        "wp_to_uniprot": {},          # WP_xxx -> uniprot_accession
        "mapping_done":  False,
        "tm_features":   {},          # uniprot_accession -> {"helices": n, "strands": n}
        "nonext_tm":     {},          # uniprot_accession -> {"helices": n, "strands": n}
        "nonext_done":   False,
    }


def save_progress(p):
    with open(PROGRESS_F, "w") as f:
        json.dump(p, f, indent=2)


def http_get(url, accept="application/json", retries=MAX_RETRIES):
    for attempt in range(1, retries + 1):
        try:
            req = urllib.request.Request(url)
            req.add_header("Accept", accept)
            with urllib.request.urlopen(req, timeout=60) as resp:
                return resp.read().decode()
        except urllib.error.HTTPError as e:
            if e.code == 429:
                wait = 60
                print(f"    Rate limited — waiting {wait}s ...")
                time.sleep(wait)
            elif e.code in (400, 404):
                return None          # not found / bad request — skip
            else:
                wait = 2 ** attempt
                print(f"    HTTP {e.code} on attempt {attempt}/{retries} — retrying in {wait}s")
                time.sleep(wait)
        except Exception as e:
            wait = 2 ** attempt
            print(f"    Error on attempt {attempt}/{retries}: {e} — retrying in {wait}s")
            time.sleep(wait)
    return None


def http_post(url, data_dict, retries=MAX_RETRIES):
    body = urllib.parse.urlencode(data_dict).encode()
    for attempt in range(1, retries + 1):
        try:
            req = urllib.request.Request(url, data=body, method="POST")
            req.add_header("Content-Type", "application/x-www-form-urlencoded")
            req.add_header("Accept", "application/json")
            with urllib.request.urlopen(req, timeout=60) as resp:
                return json.loads(resp.read().decode())
        except urllib.error.HTTPError as e:
            wait = 2 ** attempt
            print(f"    POST HTTP {e.code} attempt {attempt}/{retries} — retrying in {wait}s")
            time.sleep(wait)
        except Exception as e:
            wait = 2 ** attempt
            print(f"    POST error attempt {attempt}/{retries}: {e} — retrying in {wait}s")
            time.sleep(wait)
    return None


def submit_id_mapping(ids):
    resp = http_post(f"{BASE_UNIPROT}/idmapping/run", {
        "from": "RefSeq_Protein",
        "to":   "UniProtKB",
        "ids":  ",".join(ids),
    })
    return resp.get("jobId") if resp else None


def poll_id_mapping(job_id):
    for i in range(1, MAX_POLL + 1):
        raw = http_get(f"{BASE_UNIPROT}/idmapping/status/{job_id}")
        if raw is None:
            time.sleep(POLL_INTERVAL)
            continue
        data = json.loads(raw)
        status = data.get("jobStatus", "")
        if status == "FINISHED" or "results" in data:
            return True
        elif status in ("ERROR", "FAILED"):
            print(f"    Mapping job {job_id} failed: {status}")
            return False
        if i % 6 == 0:
            print(f"    [{job_id}] status={status} (poll {i})")
        time.sleep(POLL_INTERVAL)
    print(f"    Mapping job {job_id} timed out.")
    return False


def fetch_mapping_results(job_id):
    url  = f"{BASE_UNIPROT}/idmapping/results/{job_id}?format=json&size=500"
    raw  = http_get(url)
    if not raw:
        return {}
    data    = json.loads(raw)
    mapping = {}
    for entry in data.get("results", []):
        from_id = entry.get("from", "")
        to_acc  = entry.get("to", {})
        if isinstance(to_acc, dict):
            to_acc = to_acc.get("primaryAccession", "")
        if from_id and to_acc:
            mapping[from_id] = to_acc
    return mapping


def map_wp_ids(wp_ids, progress):
    already  = set(progress["wp_to_uniprot"].keys())
    to_map   = [x for x in wp_ids if x not in already]
    print(f"  WP_ IDs to map: {len(to_map):,} (already mapped: {len(already):,})")

    batches = [to_map[i:i+MAP_BATCH] for i in range(0, len(to_map), MAP_BATCH)]
    for bi, batch in enumerate(batches):
        print(f"\n  ID-mapping batch {bi+1}/{len(batches)} ({len(batch)} IDs) ...")
        job_id = submit_id_mapping(batch)
        if not job_id:
            print("    Failed to submit — skipping batch.")
            continue
        print(f"    Job ID: {job_id}")
        ok = poll_id_mapping(job_id)
        if not ok:
            continue
        mapping = fetch_mapping_results(job_id)
        print(f"    Mapped: {len(mapping):,}/{len(batch)}")
        progress["wp_to_uniprot"].update(mapping)
        save_progress(progress)
        time.sleep(2)

    progress["mapping_done"] = True
    save_progress(progress)
    return progress["wp_to_uniprot"]


def parse_tm_features(features_list):
    # Counts TM helices and TM strands from the UniProt features
    helices = strands = 0
    for feat in features_list:
        t = feat.get("type", "").lower()
        if t == "transmembrane":
            desc = feat.get("description", "").lower()
            if "beta" in desc or "strand" in desc:
                strands += 1
            else:
                helices += 1   # default: alpha helix
    return helices, strands


def query_tm_features_batch(accessions, progress_cache):
    to_query = [a for a in accessions if a not in progress_cache]
    if not to_query:
        return {}

    results = {}
    batches = [to_query[i:i+BATCH_SIZE] for i in range(0, len(to_query), BATCH_SIZE)]

    for bi, batch in enumerate(batches):
        query = " OR ".join(f"accession:{a}" for a in batch)
        url   = (f"{BASE_UNIPROT}/uniprotkb/search"
                 f"?query={urllib.parse.quote(query)}"
                 f"&fields=accession,ft_transmem"
                 f"&format=json&size={BATCH_SIZE}")
        raw   = http_get(url)
        if not raw:
            print(f"    Batch {bi+1}: query failed — skipping {len(batch)} accessions.")
            continue

        data = json.loads(raw)
        for entry in data.get("results", []):
            acc      = entry.get("primaryAccession", "")
            features = entry.get("features", [])
            h, s     = parse_tm_features(features)
            results[acc] = {"helices": h, "strands": s}

        # Mark not-found as queried with 0
        for a in batch:
            if a not in results:
                results[a] = {"helices": 0, "strands": 0}

        time.sleep(1)

    return results


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    progress = load_progress()

    # Step 1: extremophile sample
    print("Loading prokaryote extremophile sample ...")
    prok_records = read_fasta_ids(PROK_SAMPLE)   # (WP_id, class)
    wp_ids       = [pid for pid, _ in prok_records]
    pid2class    = {pid: cls for pid, cls in prok_records}
    print(f"  {len(wp_ids):,} sequences across {len(set(pid2class.values()))} classes")

    # Step 2: map WP_ IDs to UniProt
    if not progress["mapping_done"]:
        print("\nMapping WP_ IDs to UniProtKB ...")
        wp2uni = map_wp_ids(wp_ids, progress)
    else:
        wp2uni = progress["wp_to_uniprot"]
        print(f"\nMapping already done: {len(wp2uni):,} WP_ to UniProt mappings loaded.")

    mapped_uni = list(set(wp2uni.values()))
    print(f"  Unique UniProt accessions from mapping: {len(mapped_uni):,}")

    # Step 3: TM features of the extremophiles
    print("\nQuerying TM features for extremophile mapped accessions ...")
    cached = progress["tm_features"]
    new_tm = query_tm_features_batch(mapped_uni, cached)
    cached.update(new_tm)
    progress["tm_features"] = cached
    save_progress(progress)
    print(f"  TM features retrieved for {len(cached):,} accessions.")

    # Step 4: sample and query the non-extremophiles
    if not progress["nonext_done"]:
        print("\nSampling non-extremophile sequences ...")
        nonext_ids = read_nonextrem_sample(NONEXT_FAA)

        print("Querying TM features for non-extremophile accessions ...")
        nonext_cached = progress["nonext_tm"]
        nonext_new    = query_tm_features_batch(nonext_ids, nonext_cached)
        nonext_cached.update(nonext_new)
        progress["nonext_tm"]  = nonext_cached
        progress["nonext_done"] = True
        save_progress(progress)
    else:
        nonext_ids    = list(progress["nonext_tm"].keys())
        nonext_cached = progress["nonext_tm"]
        print(f"\nNon-extremophile TM data already done: {len(nonext_cached):,} accessions.")

    # Step 5: write the extremophile table
    prok_tsv = OUT_DIR / "prokaryote_tm.tsv"
    with open(prok_tsv, "w") as f:
        f.write("wp_id\textremophile_class\tuniprot_accession\thas_tm\tn_tm_helices\tn_tm_strands\n")
        for wp_id, cls in prok_records:
            uni  = wp2uni.get(wp_id, "")
            feat = cached.get(uni, {"helices": 0, "strands": 0}) if uni else {"helices": 0, "strands": 0}
            has  = 1 if (feat["helices"] + feat["strands"] > 0) else 0
            f.write(f"{wp_id}\t{cls}\t{uni}\t{has}\t{feat['helices']}\t{feat['strands']}\n")
    print(f"\nWrote: {prok_tsv}")

    # Step 6: write the non-extremophile table
    nonext_tsv = OUT_DIR / "nonextremophile_tm.tsv"
    with open(nonext_tsv, "w") as f:
        f.write("uniprot_accession\thas_tm\tn_tm_helices\tn_tm_strands\n")
        for acc in progress["nonext_tm"]:
            feat = progress["nonext_tm"][acc]
            has  = 1 if (feat["helices"] + feat["strands"] > 0) else 0
            f.write(f"{acc}\t{has}\t{feat['helices']}\t{feat['strands']}\n")
    print(f"Wrote: {nonext_tsv}")

    # Step 7: summary
    print("\nTM summary:")

    nonext_total  = len(progress["nonext_tm"])
    nonext_tm_n   = sum(1 for v in progress["nonext_tm"].values()
                        if v["helices"] + v["strands"] > 0)
    nonext_frac   = nonext_tm_n / nonext_total if nonext_total else 0.0
    print(f"  Non-extremophile: {nonext_tm_n}/{nonext_total} TM = {nonext_frac:.3f}")

    # Per-class extremophile (only proteins that mapped to UniProt)
    class_stats = {}
    for cls in CLASS_NAMES:
        cls_ids     = [pid for pid, c in prok_records if c == cls]
        n_total     = len(cls_ids)
        n_mapped    = sum(1 for pid in cls_ids if pid in wp2uni)
        n_tm        = sum(1 for pid in cls_ids
                          if pid in wp2uni
                          and cached.get(wp2uni[pid], {"helices":0,"strands":0})["helices"]
                            + cached.get(wp2uni[pid], {"helices":0,"strands":0})["strands"] > 0)
        frac_mapped = n_tm / n_mapped if n_mapped else 0.0
        class_stats[cls] = {
            "n_sampled": n_total,
            "n_mapped_to_uniprot": n_mapped,
            "n_with_tm": n_tm,
            "tm_fraction_of_mapped": frac_mapped,
        }
        print(f"  {cls:15s}: {n_tm}/{n_mapped} mapped TM = {frac_mapped:.3f}  "
              f"({n_mapped}/{n_total} IDs mapped)")

    summary = {
        "non_extremophile": {
            "n_queried": nonext_total,
            "n_with_tm": nonext_tm_n,
            "tm_fraction": nonext_frac,
        },
        "extremophile_classes": class_stats,
    }
    with open(OUT_DIR / "tm_summary.json", "w") as f:
        json.dump(summary, f, indent=2)
    print(f"\nSaved summary to {OUT_DIR / 'tm_summary.json'}")
    print("\nDone.")


if __name__ == "__main__":
    main()
