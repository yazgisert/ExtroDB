#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import math
import re
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple


BASE = Path("/Users/yazgisert/Documents/Playground/thesis_data_collection")
EXISTING_SCREEN_LOG = BASE / "paper_screening_log_tempura_bacdive.tsv"
DEFAULT_OUTPUT = BASE / "paper_screening_log_tempura_broad_second_pass.tsv"

THERMOPHILE_MIN = 45.0
HYPERTHERMOPHILE_MIN = 80.0


def normalize_text(value: str) -> str:
    return re.sub(r"[^a-z0-9]", "", value.lower())


def detect_delimiter(path: Path) -> str:
    with path.open("r", encoding="utf-8", newline="") as handle:
        sample = handle.read(2048)
    if sample.count("\t") >= sample.count(","):
        return "\t"
    return ","


def parse_float(value: str) -> Optional[float]:
    if value is None:
        return None
    text = str(value).strip()
    if not text:
        return None
    match = re.search(r"-?\d+(?:\.\d+)?", text)
    if not match:
        return None
    try:
        return float(match.group(0))
    except ValueError:
        return None


def read_table(path: Path) -> List[Dict[str, str]]:
    delimiter = detect_delimiter(path)
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter=delimiter))


def find_column(fieldnames: Sequence[str], candidates: Sequence[str]) -> str:
    normalized = {name.lower().strip(): name for name in fieldnames}
    for candidate in candidates:
        exact = normalized.get(candidate.lower().strip())
        if exact:
            return exact
    for name in fieldnames:
        low = name.lower().strip()
        for candidate in candidates:
            if candidate.lower().strip() in low:
                return name
    raise KeyError(f"Could not find any of these columns: {candidates}")


def existing_screened_keys(rows: Iterable[Dict[str, str]]) -> set[Tuple[str, str]]:
    keys = set()
    for row in rows:
        taxon = normalize_text(row.get("taxon_name_candidate", ""))
        strain = normalize_text(row.get("strain_candidate", ""))
        if taxon and strain:
            keys.add((taxon, strain))
    return keys


def build_output_row(
    taxon_name: str,
    strain: str,
    topt: float,
) -> Dict[str, str]:
    label = "hyperthermophile" if topt >= HYPERTHERMOPHILE_MIN else "thermophile"
    return {
        "source_type": "TEMPURA",
        "citation_or_resource_url": "https://togodb.org/db/tempura",
        "focus_category": "thermophile",
        "taxon_name_candidate": taxon_name,
        "strain_candidate": strain,
        "key_evidence": f"TEMPURA Topt_average = {topt:g} C",
        "extracted_ogt_celsius": f"{topt:g}",
        "extracted_nacl_optimum_M": "",
        "extracted_label": label,
        "evidence_level_candidate": "curated_database",
        "follow_up_action": "validate_in_BacDive",
        "notes": "Broad TEMPURA second-pass candidate; not yet BacDive-validated.",
    }


def score_row(taxon_name: str, strain: str, topt: float) -> Tuple[int, float, str, str]:
    # Prefer stable-looking strain-aware rows first, then hottest candidates.
    has_collection_style = bool(re.match(r"^(DSM|JCM|NBRC|ATCC|IFO|KCTC|CGMCC|LMG|CIP)\s*\d+", strain, re.I))
    has_clean_strain = 0 if strain.strip() else 1
    collection_penalty = 0 if has_collection_style else 1
    return (has_clean_strain, collection_penalty, -topt, taxon_name.lower())


def prepare_candidates(
    tempura_rows: List[Dict[str, str]],
    screened_keys: set[Tuple[str, str]],
    species_col: str,
    strain_col: str,
    topt_col: str,
    limit: int,
) -> List[Dict[str, str]]:
    prepared: List[Tuple[Tuple[int, float, str, str], Dict[str, str]]] = []
    seen = set(screened_keys)

    for row in tempura_rows:
        taxon_name = (row.get(species_col) or "").strip()
        strain = (row.get(strain_col) or "").strip()
        topt = parse_float(row.get(topt_col, ""))

        if not taxon_name or not strain or topt is None:
            continue
        if topt < THERMOPHILE_MIN:
            continue

        key = (normalize_text(taxon_name), normalize_text(strain))
        if not key[0] or not key[1]:
            continue
        if key in seen:
            continue

        out_row = build_output_row(taxon_name, strain, topt)
        prepared.append((score_row(taxon_name, strain, topt), out_row))
        seen.add(key)

    prepared.sort(key=lambda item: item[0])
    rows = [row for _, row in prepared]
    if limit > 0:
        rows = rows[:limit]
    return rows


def write_tsv(path: Path, rows: List[Dict[str, str]]) -> None:
    if not rows:
        raise SystemExit("No new TEMPURA second-pass candidates were found.")
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Prepare a broader TEMPURA thermophile second-pass screening log."
    )
    parser.add_argument("tempura_export", help="Path to the broader TEMPURA export CSV/TSV.")
    parser.add_argument(
        "--output",
        default=str(DEFAULT_OUTPUT),
        help="Output TSV path for the broader second-pass screening log.",
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=50,
        help="Maximum number of new candidates to emit. Use 0 for no limit.",
    )
    args = parser.parse_args()

    export_path = Path(args.tempura_export)
    output_path = Path(args.output)

    tempura_rows = read_table(export_path)
    if not tempura_rows:
        raise SystemExit("TEMPURA export has no rows.")

    existing_rows = read_table(EXISTING_SCREEN_LOG)
    screened_keys = existing_screened_keys(existing_rows)

    fieldnames = list(tempura_rows[0].keys())
    species_col = find_column(fieldnames, ["Genus_and_species", "genus_and_species", "species"])
    strain_col = find_column(fieldnames, ["Strain", "strain"])
    topt_col = find_column(
        fieldnames,
        ["Topt_average(°C)", "Topt_average", "topt_average", "Topt_ave", "topt_ave"],
    )

    output_rows = prepare_candidates(
        tempura_rows=tempura_rows,
        screened_keys=screened_keys,
        species_col=species_col,
        strain_col=strain_col,
        topt_col=topt_col,
        limit=args.limit,
    )
    write_tsv(output_path, output_rows)

    hyper = sum(1 for row in output_rows if row["extracted_label"] == "hyperthermophile")
    thermo = sum(1 for row in output_rows if row["extracted_label"] == "thermophile")
    print(f"output={output_path}")
    print(f"rows={len(output_rows)}")
    print(f"thermophile={thermo}")
    print(f"hyperthermophile={hyper}")


if __name__ == "__main__":
    main()
