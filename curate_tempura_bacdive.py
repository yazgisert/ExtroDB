#!/usr/bin/env python3

from __future__ import annotations

import json
import math
import re
import sys
import time
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Tuple
from urllib.error import HTTPError, URLError
from urllib.parse import quote

from bacdive_client import (
    BACDIVE_API_BASE,
    BACDIVE_USER_AGENT,
    ensure_list,
    fetch_json as request_json,
    portal_url,
)
from common import read_tsv, write_tsv

API_BASE = BACDIVE_API_BASE
USER_AGENT = BACDIVE_USER_AGENT
REQUEST_DELAY_SECONDS = 0.05
CORE_FIELDNAMES = [
    "taxon_name",
    "ncbi_taxid",
    "strain",
    "thermophile",
    "hyperthermophile",
    "halophile",
    "halophile_class",
    "ogt_celsius",
    "nacl_optimum_M",
    "evidence_type",
    "evidence_source",
    "include_in_core",
    "notes",
]


def normalize_text(value: str) -> str:
    # Strain names vary a lot in punctuation and spacing, so we normalize them
    # before matching.
    return re.sub(r"[^a-z0-9]", "", value.lower())


def is_collection_number(value: str) -> bool:
    return bool(
        re.match(r"^(DSM|JCM|NBRC|ATCC|CIP|CCUG|NCIMB|CGMCC|IFO|KCTC|LMG)\s*\d+", value.strip(), re.I)
    )


def try_fetch_json(url: str) -> Optional[Dict[str, Any]]:
    # Missing BacDive endpoints are treated as "no result" rather than hard
    # failures so one bad lookup does not stop the whole curation pass.
    try:
        time.sleep(REQUEST_DELAY_SECONDS)
        return request_json(url, user_agent=USER_AGENT)
    except HTTPError as exc:
        if exc.code == 404:
            return None
        raise
    except URLError:
        return None


def parse_taxon_name(full_name: str) -> Tuple[str, str]:
    parts = full_name.strip().split()
    if len(parts) < 2:
        raise ValueError(f"Expected species-level taxon name, got: {full_name!r}")
    return parts[0], parts[1]


def parse_numeric_temperature(value: Any) -> Optional[float]:
    if value is None:
        return None
    if isinstance(value, (int, float)):
        if math.isnan(value):
            return None
        return float(value)
    text = str(value).strip()
    match = re.search(r"-?\d+(?:\.\d+)?", text)
    if not match:
        return None
    return float(match.group(0))


def extract_ncbi_taxid(record: Dict[str, Any]) -> Optional[str]:
    general = record.get("General", {})
    ncbi_field = general.get("NCBI tax id")
    entries = ensure_list(ncbi_field)
    best = None
    for entry in entries:
        if isinstance(entry, dict):
            taxid = entry.get("NCBI tax id")
            level = str(entry.get("Matching level", "")).lower()
        else:
            taxid = entry
            level = ""
        if taxid is None:
            continue
        if "species" in level:
            return str(taxid)
        if best is None:
            best = str(taxid)
    return best


def extract_strain_variants(record: Dict[str, Any]) -> List[str]:
    # We collect strain aliases from multiple BacDive sections because the
    # candidate table and BacDive do not always name the same strain the same way.
    variants: List[str] = []

    name_section = record.get("Name and taxonomic classification", {})
    strain_designation = name_section.get("strain designation")
    for value in ensure_list(strain_designation):
        if isinstance(value, dict):
            raw = value.get("strain designation")
        else:
            raw = value
        if raw:
            parts = re.split(r"[,;]", str(raw))
            variants.extend(part.strip() for part in parts if part.strip())

    literature = record.get("Literature", {})
    culture_numbers = literature.get("culture collection no.")
    if culture_numbers:
        parts = re.split(r"[,;]", str(culture_numbers))
        variants.extend(part.strip() for part in parts if part.strip())

    general = record.get("General", {})
    dsm_number = general.get("DSM-Number")
    if dsm_number:
        variants.append(f"DSM {dsm_number}")

    for history_item in ensure_list(general.get("strain history")):
        if isinstance(history_item, dict):
            history = history_item.get("history")
            if history:
                variants.append(str(history))

    seen = set()
    unique: List[str] = []
    for variant in variants:
        if variant not in seen:
            unique.append(variant)
            seen.add(variant)
    return unique


def species_exact_match(record: Dict[str, Any], candidate_taxon: str) -> bool:
    name_section = record.get("Name and taxonomic classification", {})
    species = str(name_section.get("species", "")).strip().lower()
    return species == candidate_taxon.strip().lower()


def score_record(record: Dict[str, Any], candidate_taxon: str, candidate_strain: str) -> Tuple[int, str]:
    # Matching is intentionally transparent: higher scores mean better species
    # and strain agreement, and the reason string shows why a record won.
    score = 0
    reasons: List[str] = []
    candidate_norm = normalize_text(candidate_strain)

    if species_exact_match(record, candidate_taxon):
        score += 20
        reasons.append("species")

    variants = extract_strain_variants(record)
    variant_norms = [(variant, normalize_text(variant)) for variant in variants]

    exact_match = False
    partial_match = False

    for variant, norm in variant_norms:
        if not norm:
            continue
        if candidate_norm and norm == candidate_norm:
            score += 100
            reasons.append(f"exact:{variant}")
            exact_match = True
            break

    if not exact_match and candidate_norm:
        for variant, norm in variant_norms:
            if not norm:
                continue
            if candidate_norm in norm or norm in candidate_norm:
                score += 75
                reasons.append(f"partial:{variant}")
                partial_match = True
                break

    if not exact_match and not partial_match and candidate_norm:
        for variant in variants:
            if candidate_strain.lower() in variant.lower():
                score += 60
                reasons.append(f"text:{variant}")
                partial_match = True
                break

    if name_section_type_strain(record):
        score += 3
        reasons.append("type")

    return score, ",".join(reasons)


def name_section_type_strain(record: Dict[str, Any]) -> bool:
    name_section = record.get("Name and taxonomic classification", {})
    value = str(name_section.get("type strain", "")).strip().lower()
    return value == "yes"


def extract_temperature_entries(record: Dict[str, Any]) -> List[Dict[str, Any]]:
    section = record.get("Culture and growth conditions", {})
    temp_field = section.get("culture temp")
    entries: List[Dict[str, Any]] = []
    for item in ensure_list(temp_field):
        if isinstance(item, dict):
            entries.append(item)
    return entries


def summarize_temperature(record: Dict[str, Any]) -> Dict[str, Any]:
    # BacDive can contain optimum, growth, minimum, or maximum temperatures.
    # We prefer optimum values first and fall back only when needed.
    entries = extract_temperature_entries(record)
    optimum_values: List[float] = []
    growth_values: List[float] = []
    max_values: List[float] = []
    min_values: List[float] = []

    for entry in entries:
        temp = parse_numeric_temperature(entry.get("temperature"))
        if temp is None:
            continue
        entry_type = str(entry.get("type", "")).strip().lower()
        if "opt" in entry_type:
            optimum_values.append(temp)
        elif "max" in entry_type:
            max_values.append(temp)
        elif "min" in entry_type:
            min_values.append(temp)
        else:
            growth_values.append(temp)

    chosen = None
    chosen_type = None
    if optimum_values:
        chosen = optimum_values[0]
        chosen_type = "optimum"
    elif growth_values:
        chosen = max(growth_values)
        chosen_type = "growth"
    elif max_values:
        chosen = max(max_values)
        chosen_type = "maximum"
    elif min_values:
        chosen = max(min_values)
        chosen_type = "minimum"

    return {
        "entries": entries,
        "chosen_temperature": chosen,
        "chosen_type": chosen_type,
        "optimum_values": optimum_values,
        "growth_values": growth_values,
        "maximum_values": max_values,
        "minimum_values": min_values,
    }


def compare_tempura_to_bacdive(tempura_ogt: Optional[float], bacdive_temp: Optional[float]) -> str:
    if tempura_ogt is None or bacdive_temp is None:
        return "unknown"
    if abs(tempura_ogt - bacdive_temp) <= 3:
        return "agrees"
    if tempura_ogt >= 80 and bacdive_temp >= 80:
        return "supports_same_class"
    if tempura_ogt >= 45 and bacdive_temp >= 45:
        return "supports_same_class"
    return "differs"


def make_note(
    taxon: str,
    strain: str,
    match_status: str,
    bacdive_id: Optional[int],
    temperature_summary: Dict[str, Any],
    tempura_ogt: Optional[float],
    agreement: str,
) -> str:
    if match_status == "no_bacdive_match":
        return "No BacDive strain entry found for this TEMPURA candidate."

    chosen = temperature_summary.get("chosen_temperature")
    chosen_type = temperature_summary.get("chosen_type")
    record_url = portal_url(bacdive_id) if bacdive_id is not None else ""

    pieces = [f"BacDive match for {taxon} {strain}".strip()]
    if bacdive_id is not None:
        pieces.append(f"BacDive ID {bacdive_id}")
    if chosen is not None and chosen_type:
        pieces.append(f"{chosen_type} temperature evidence {format_temperature(chosen)} C")
    else:
        pieces.append("no usable culture_temp value found")
    if agreement == "agrees":
        pieces.append("BacDive agrees closely with TEMPURA")
    elif agreement == "supports_same_class":
        pieces.append("BacDive supports the same thermophile class as TEMPURA")
    elif agreement == "differs":
        pieces.append("BacDive temperature differs from TEMPURA and should be reviewed")
    else:
        pieces.append("TEMPURA to BacDive comparison could not be determined automatically")
    if record_url:
        pieces.append(record_url)
    return ". ".join(pieces) + "."


def format_temperature(value: Optional[float]) -> str:
    if value is None:
        return ""
    if float(value).is_integer():
        return str(int(value))
    return f"{value:.1f}".rstrip("0").rstrip(".")


def classify_include(
    exact_match: bool,
    temperature_summary: Dict[str, Any],
    tempura_ogt: Optional[float],
) -> str:
    # A row enters the core only if the BacDive match is strong enough and the
    # temperature evidence supports thermophily.
    chosen = temperature_summary.get("chosen_temperature")
    chosen_type = temperature_summary.get("chosen_type")
    agreement = compare_tempura_to_bacdive(tempura_ogt, chosen)

    if chosen is None:
        return "review"
    if chosen < 45:
        return "no"
    if not exact_match:
        return "review"
    if agreement == "differs":
        return "review"
    return "yes"


def classify_labels(chosen_temperature: Optional[float]) -> Tuple[str, str]:
    if chosen_temperature is None or chosen_temperature < 45:
        return "no", "no"
    if chosen_temperature >= 80:
        return "yes", "yes"
    return "yes", "no"


def build_bacdive_queries(taxon: str, strain: str) -> List[str]:
    genus, species = parse_taxon_name(taxon)
    queries: List[str] = []
    if is_collection_number(strain):
        queries.append(f"{API_BASE}/culturecollectionno/{quote(strain)}")
    queries.append(f"{API_BASE}/taxon/{quote(genus)}/{quote(species)}")
    return queries


def fetch_candidate_ids(taxon: str, strain: str) -> List[int]:
    ids: List[int] = []
    seen = set()
    for url in build_bacdive_queries(taxon, strain):
        payload = try_fetch_json(url)
        if not payload:
            continue
        for value in payload.get("results", []):
            try:
                bacdive_id = int(value)
            except (TypeError, ValueError):
                continue
            if bacdive_id not in seen:
                seen.add(bacdive_id)
                ids.append(bacdive_id)
    return ids


def fetch_record(bacdive_id: int) -> Optional[Dict[str, Any]]:
    payload = try_fetch_json(f"{API_BASE}/fetch/{bacdive_id}")
    if not payload:
        return None
    results = payload.get("results", {})
    record = results.get(str(bacdive_id))
    if not isinstance(record, dict):
        return None
    return record


def iter_scored_records(
    taxon: str,
    strain: str,
) -> Iterable[Tuple[int, Dict[str, Any], int, str]]:
    # Separating iteration from selection keeps the "score every candidate"
    # logic easy to read and test.
    for bacdive_id in fetch_candidate_ids(taxon, strain):
        record = fetch_record(bacdive_id)
        if not record:
            continue
        score, reason = score_record(record, taxon, strain)
        yield bacdive_id, record, score, reason


def choose_best_record(taxon: str, strain: str) -> Tuple[Optional[int], Optional[Dict[str, Any]], int, str]:
    # We keep only the highest-scoring BacDive record for each TEMPURA row.
    best_match = max(
        iter_scored_records(taxon, strain),
        key=lambda item: item[2],
        default=(None, None, -1, ""),
    )
    return best_match


def should_accept_match(score: int) -> bool:
    return score >= 90


def curate_row(row: Dict[str, str]) -> Dict[str, Any]:
    # This is the main per-row workflow:
    # 1. find the best BacDive record,
    # 2. summarize temperature evidence,
    # 3. decide inclusion,
    # 4. build screening/core output rows.
    taxon = row["taxon_name_candidate"].strip()
    strain = row["strain_candidate"].strip()
    tempura_ogt = parse_numeric_temperature(row.get("extracted_ogt_celsius"))

    bacdive_id, record, score, score_reason = choose_best_record(taxon, strain)
    if not record or bacdive_id is None or not should_accept_match(score):
        row_copy = dict(row)
        row_copy["notes"] = "No BacDive strain entry found for this candidate; excluded from first BacDive-backed core table."
        return {
            "screening_row": row_copy,
            "core_row": None,
            "status": "no_bacdive_match",
            "taxon": taxon,
            "strain": strain,
        }

    temperature_summary = summarize_temperature(record)
    chosen_temp = temperature_summary["chosen_temperature"]
    chosen_type = temperature_summary["chosen_type"]
    agreement = compare_tempura_to_bacdive(tempura_ogt, chosen_temp)
    exact_match = score >= 120
    include_in_core = classify_include(exact_match, temperature_summary, tempura_ogt)
    thermophile, hyperthermophile = classify_labels(chosen_temp)
    note = make_note(
        taxon=taxon,
        strain=strain,
        match_status="bacdive_match",
        bacdive_id=bacdive_id,
        temperature_summary=temperature_summary,
        tempura_ogt=tempura_ogt,
        agreement=agreement,
    )

    screening_row = dict(row)
    screening_row["notes"] = note

    taxon_name = str(record.get("Name and taxonomic classification", {}).get("species", taxon))
    core_row = {
        "taxon_name": taxon_name,
        "ncbi_taxid": extract_ncbi_taxid(record) or "",
        "strain": strain,
        "thermophile": thermophile,
        "hyperthermophile": hyperthermophile,
        "halophile": "no",
        "halophile_class": "",
        "ogt_celsius": format_temperature(chosen_temp) if chosen_temp is not None else "",
        "nacl_optimum_M": "",
        "evidence_type": "curated_database",
        "evidence_source": portal_url(bacdive_id),
        "include_in_core": include_in_core,
        "notes": (
            f"Matched BacDive ID {bacdive_id} ({score_reason}). "
            f"Used {chosen_type or 'no'} temperature evidence from BacDive."
        ),
    }

    return {
        "screening_row": screening_row,
        "core_row": core_row,
        "status": include_in_core,
        "taxon": taxon,
        "strain": strain,
    }


def read_input_rows(path: Path) -> List[Dict[str, str]]:
    return read_tsv(path)


def main() -> int:
    if len(sys.argv) != 4:
        print(
            "Usage: curate_tempura_bacdive.py INPUT_TEMPLATE_TSV OUTPUT_SCREENING_TSV OUTPUT_CORE_TSV",
            file=sys.stderr,
        )
        return 1

    input_path = Path(sys.argv[1])
    screening_path = Path(sys.argv[2])
    core_path = Path(sys.argv[3])

    input_rows = read_input_rows(input_path)

    curated = [curate_row(row) for row in input_rows]
    screening_rows = [entry["screening_row"] for entry in curated]
    core_rows = [entry["core_row"] for entry in curated if entry["core_row"] is not None]

    write_tsv(screening_path, screening_rows, screening_rows[0].keys())
    if core_rows:
        write_tsv(core_path, core_rows, core_rows[0].keys())
    else:
        write_tsv(core_path, [], CORE_FIELDNAMES)

    excluded = [
        f"{entry['taxon']} {entry['strain']}".strip()
        for entry in curated
        if entry["core_row"] is None
    ]
    summary = {
        "input_rows": len(input_rows),
        "screening_rows": len(screening_rows),
        "core_rows": len(core_rows),
        "excluded_no_bacdive": excluded,
    }
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
