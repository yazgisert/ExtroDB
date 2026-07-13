#!/usr/bin/env python3

from __future__ import annotations

import json
import re
import ssl
import time
from typing import Any

from bacdive_client import (
    BACDIVE_API_BASE,
    BACDIVE_USER_AGENT,
    ensure_list,
    fetch_json,
    portal_url,
)
from common import BASE_DIR, read_tsv, write_tsv


BASE = BASE_DIR
CANDIDATES = BASE / "bacdive_halophile_candidates.tsv"
SCREENING_OUT = BASE / "paper_screening_log_halophile_bacdive.tsv"
CORE_OUT = BASE / "core_taxon_table_halophile_bacdive.tsv"
SHORTLIST_OUT = BASE / "bacdive_halophile_shortlist.tsv"

API_BASE = f"{BACDIVE_API_BASE}/fetch"
USER_AGENT = BACDIVE_USER_AGENT
REQUEST_DELAY = 0.08
TIMEOUT = 30
CLASS_QUOTAS = {
    "extreme_halophile": 11,
    "moderate_halophile": 11,
    "slight_halophile": 11,
}

SLIGHT_MIN = 0.2
MODERATE_MIN = 0.85
EXTREME_MIN = 3.4
BOUNDARY_EPS = 0.05

EXTREME_GENUS_RE = re.compile(
    r"^(Halo|Natr|Natro|Halalkali|Haladapt|Halarch|Halobaculum|Halobacter|Halococcus|"
    r"Haloferax|Haloarcula|Halogeometricum|Halomarina|Halomicrobium|Halopelagius|"
    r"Haloplanus|Halopiger|Halorubrum|Halostagnicola|Haloterrigena|Salinibacter)",
    re.I,
)
MODERATE_GENUS_RE = re.compile(
    r"^(Halomonas|Salinivibrio|Salinicola|Salinicoccus|Chromohalobacter|Halobacillus|"
    r"Lentibacillus|Pontibacillus|Marinococcus|Virgibacillus|Salimicrobium|"
    r"Tetragenococcus|Gracilibacillus|Thiohalobacter|Aliifodinibius|Idiomarina)",
    re.I,
)


def read_candidates() -> list[dict[str, str]]:
    return read_tsv(CANDIDATES)


def request_json(url: str) -> dict[str, Any]:
    return fetch_json(url, user_agent=USER_AGENT, timeout=TIMEOUT)


def fetch_record(bacdive_id: str) -> dict[str, Any]:
    # A short pause keeps the BacDive API requests polite and predictable.
    time.sleep(REQUEST_DELAY)
    payload = request_json(f"{API_BASE}/{bacdive_id}")
    return payload["results"][str(bacdive_id)]


def parse_numbers(text: str) -> list[float]:
    return [float(x) for x in re.findall(r"\d+(?:\.\d+)?", text)]


def parse_concentration_to_molar(text: str) -> dict[str, float | str] | None:
    # BacDive stores NaCl optima in mixed units; this converts them to a common
    # molar scale so the class thresholds are comparable.
    raw = text.strip()
    low = raw.lower().replace("–", "-")
    nums = parse_numbers(low)
    if not nums:
        return None

    if "%" in low:
        values = [(n * 10.0) / 58.44 for n in nums]
        unit = "%(w/v)"
    elif re.search(r"(?<![a-z])m(?![a-z])", low):
        values = nums
        unit = "M"
    elif "g/l" in low or "g l-1" in low or "g l−1" in low:
        values = [n / 58.44 for n in nums]
        unit = "g/L"
    else:
        return None

    if len(values) == 1:
        midpoint = values[0]
        min_v = max_v = values[0]
    else:
        min_v = min(values)
        max_v = max(values)
        midpoint = (min_v + max_v) / 2.0

    return {
        "raw": raw,
        "unit": unit,
        "midpoint_M": midpoint,
        "min_M": min_v,
        "max_M": max_v,
    }


def classify_halophile(
    conc: dict[str, float | str] | None,
) -> tuple[str | None, str, str]:
    # We classify only from the NaCl optimum. Borderline values stay in review
    # instead of being forced into a class automatically.
    if conc is None:
        return None, "review", "unparseable concentration"

    min_v = float(conc["min_M"])
    max_v = float(conc["max_M"])
    mid = float(conc["midpoint_M"])

    if mid <= SLIGHT_MIN:
        return None, "no", "below halophile threshold"

    if max_v < MODERATE_MIN:
        label = "slight_halophile"
    elif min_v >= MODERATE_MIN and max_v < EXTREME_MIN:
        label = "moderate_halophile"
    elif min_v >= EXTREME_MIN:
        label = "extreme_halophile"
    else:
        return None, "review", "range crosses class boundary"

    for boundary in (SLIGHT_MIN, MODERATE_MIN, EXTREME_MIN):
        if abs(mid - boundary) <= BOUNDARY_EPS:
            return label, "review", "near class threshold"

    return label, "yes", "clear optimum"


def label_to_class(label: str | None) -> str:
    return "" if label is None else label.replace("_halophile", "")


def extract_ncbi_taxid(record: dict[str, Any]) -> str:
    field = record.get("General", {}).get("NCBI tax id")
    for entry in ensure_list(field):
        if isinstance(entry, dict) and entry.get("Matching level", "").lower() == "species":
            return str(entry.get("NCBI tax id", ""))
        if not isinstance(entry, dict):
            return str(entry)
    for entry in ensure_list(field):
        if isinstance(entry, dict):
            val = entry.get("NCBI tax id")
            if val is not None:
                return str(val)
    return ""


def best_strain(candidate_row: dict[str, str], record: dict[str, Any]) -> str:
    designation = candidate_row["designation_candidate"].strip()
    if designation:
        return designation
    name_section = record.get("Name and taxonomic classification", {})
    strain = name_section.get("strain designation")
    if strain:
        return str(strain).strip()
    numbers = candidate_row["strain_number_candidate"].strip()
    if numbers:
        return numbers.split(",")[0].strip()
    return ""


def extract_optimum_halophily(record: dict[str, Any]) -> list[dict[str, Any]]:
    halophily = record.get("Physiology and metabolism", {}).get("halophily")
    return [
        item
        for item in ensure_list(halophily)
        if isinstance(item, dict)
        and str(item.get("salt", "")).lower() == "nacl"
        and str(item.get("tested relation", "")).lower() == "optimum"
        and item.get("concentration")
    ]


def choose_entry_with_simplest_concentration(
    entries: list[dict[str, Any]],
) -> dict[str, Any] | None:
    return min(
        entries,
        key=lambda entry: len(parse_numbers(str(entry.get("concentration", "")))),
        default=None,
    )


def classify_record(record: dict[str, Any]) -> dict[str, Any]:
    # This bundles the main interpretation of one BacDive record so selection
    # code can reuse it without repeating the same steps.
    optimum_entries = extract_optimum_halophily(record)
    chosen_entry = choose_entry_with_simplest_concentration(optimum_entries)
    concentration = parse_concentration_to_molar(chosen_entry["concentration"]) if chosen_entry else None
    label, include, reason = classify_halophile(concentration)
    return {
        "record": record,
        "chosen_entry": chosen_entry,
        "conc": concentration,
        "label": label,
        "include": include,
        "reason": reason,
    }


def species_pool_key(row: dict[str, str]) -> str:
    species = row["species"]
    return species.split()[0] if species else ""


def shortlist_rows(rows: list[dict[str, str]]) -> dict[str, list[dict[str, str]]]:
    # The pools bias the first-pass core toward obvious halophile genera while
    # still leaving a general fallback bucket.
    type_rows = [row for row in rows if row["is_type_strain"] == "yes"]
    extreme = [row for row in type_rows if EXTREME_GENUS_RE.search(species_pool_key(row))]
    moderate = [row for row in type_rows if MODERATE_GENUS_RE.search(species_pool_key(row))]
    used = {r["bacdive_id"] for r in extreme + moderate}
    general = [row for row in type_rows if row["bacdive_id"] not in used]
    return {
        "extreme": extreme,
        "moderate": moderate,
        "general": general,
        "all": type_rows,
    }


def select_balanced_core(
    rows: list[dict[str, str]],
) -> tuple[list[tuple[dict[str, str], dict[str, Any]]], dict[str, int]]:
    # The goal is a balanced core: try to fill each halophile class quota once,
    # without selecting the same BacDive ID twice.
    pools = shortlist_rows(rows)
    cache: dict[str, dict[str, Any]] = {}
    quotas = dict(CLASS_QUOTAS)
    selected: list[tuple[dict[str, str], dict[str, Any]]] = []
    selected_ids: set[str] = set()

    target_orders = {
        "extreme_halophile": pools["extreme"] + pools["moderate"] + pools["general"] + pools["all"],
        "moderate_halophile": pools["moderate"] + pools["general"] + pools["extreme"] + pools["all"],
        "slight_halophile": pools["general"] + pools["moderate"] + pools["extreme"] + pools["all"],
    }

    for target_label, ordered_rows in target_orders.items():
        for row in ordered_rows:
            if quotas[target_label] <= 0:
                break
            bid = row["bacdive_id"]
            if bid in selected_ids:
                continue
            if bid not in cache:
                cache[bid] = classify_record(fetch_record(bid))
            data = cache[bid]
            if data["label"] != target_label:
                continue
            if data["include"] != "yes":
                continue
            selected.append((row, data))
            selected_ids.add(bid)
            quotas[target_label] -= 1

    return selected, quotas


def screening_row(candidate_row: dict[str, str], data: dict[str, Any]) -> dict[str, str]:
    record = data["record"]
    conc = data["conc"]
    bid = candidate_row["bacdive_id"]
    species = candidate_row["species"]
    strain = best_strain(candidate_row, record)
    raw = data["chosen_entry"]["concentration"] if data["chosen_entry"] else ""
    mol = conc["midpoint_M"] if conc else ""
    label = data["label"] or ""
    notes = []
    if candidate_row["is_type_strain"] == "yes":
        notes.append("type strain")
    if data["reason"]:
        notes.append(data["reason"])
    return {
        "source_type": "BacDive",
        "citation_or_resource_url": portal_url(bid),
        "focus_category": "halophile",
        "taxon_name_candidate": species,
        "strain_candidate": strain,
        "key_evidence": f"BacDive NaCl optimum = {raw}",
        "extracted_ogt_celsius": "",
        "extracted_nacl_optimum_M": f"{mol:.3f}",
        "extracted_label": label,
        "evidence_level_candidate": "curated_database",
        "follow_up_action": "transfer_to_core",
        "notes": "; ".join(notes),
    }


def core_row(candidate_row: dict[str, str], data: dict[str, Any]) -> dict[str, str]:
    record = data["record"]
    label = data["label"]
    halophile_class = label_to_class(label)
    bid = candidate_row["bacdive_id"]
    conc = data["conc"]
    return {
        "taxon_name": candidate_row["species"],
        "ncbi_taxid": extract_ncbi_taxid(record),
        "strain": best_strain(candidate_row, record),
        "thermophile": "no",
        "hyperthermophile": "no",
        "halophile": "yes",
        "halophile_class": halophile_class,
        "ogt_celsius": "",
        "nacl_optimum_M": f"{conc['midpoint_M']:.3f}",
        "evidence_type": "curated_database",
        "evidence_source": portal_url(bid),
        "include_in_core": "yes",
        "notes": (
            f"BacDive NaCl optimum {data['chosen_entry']['concentration']}; "
            f"classified as {halophile_class}; type strain {candidate_row['is_type_strain']}."
        ),
    }


def shortlist_row(candidate_row: dict[str, str], data: dict[str, Any]) -> dict[str, str]:
    record = data["record"]
    conc = data["conc"]
    raw = data["chosen_entry"]["concentration"] if data["chosen_entry"] else ""
    return {
        "bacdive_id": candidate_row["bacdive_id"],
        "species": candidate_row["species"],
        "strain": best_strain(candidate_row, record),
        "ncbi_taxid": extract_ncbi_taxid(record),
        "is_type_strain": candidate_row["is_type_strain"],
        "raw_optimum": raw,
        "nacl_optimum_M": f"{conc['midpoint_M']:.3f}" if conc else "",
        "halophile_label": data["label"] or "",
        "include_in_core": data["include"],
        "reason": data["reason"],
        "bacdive_url": portal_url(candidate_row["bacdive_id"]),
    }


def main() -> None:
    ssl._create_default_https_context = ssl._create_unverified_context
    rows = read_candidates()
    selected, quotas = select_balanced_core(rows)

    screening_rows = [screening_row(candidate, data) for candidate, data in selected]
    core_rows = [core_row(candidate, data) for candidate, data in selected]
    shortlist_rows_out = [shortlist_row(candidate, data) for candidate, data in selected]

    write_tsv(SCREENING_OUT, screening_rows, screening_rows[0].keys())
    write_tsv(CORE_OUT, core_rows, core_rows[0].keys())
    write_tsv(SHORTLIST_OUT, shortlist_rows_out, shortlist_rows_out[0].keys())

    counts = {label: sum(row["halophile_class"] == label for row in core_rows) for label in ("slight", "moderate", "extreme")}
    print(
        json.dumps(
            {
                "selected_total": len(selected),
                "class_counts": counts,
                "remaining_quotas": quotas,
                "screening_out": str(SCREENING_OUT),
                "core_out": str(CORE_OUT),
                "shortlist_out": str(SHORTLIST_OUT),
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
