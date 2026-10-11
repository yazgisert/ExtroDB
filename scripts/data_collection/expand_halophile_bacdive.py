# Second-pass halophile expansion from BacDive: takes candidates not yet in the shortlist (type strains only by
# default; known halophile genera first), classifies them with the rules of curate_halophile_bacdive.py and adds up
# to 5 per class. Output files get the suffix --batch-name (default expansion_round2).

from __future__ import annotations

import argparse
import json
import ssl
import sys
from collections import Counter
from pathlib import Path
from typing import Dict, List, Tuple

import curate_halophile_bacdive as base

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from common import BASE_DIR, read_tsv, write_tsv


BASE = BASE_DIR
CURRENT_SHORTLIST = BASE / "curation" / "source_tables" / "bacdive_halophile_shortlist.tsv"
TARGET_LABELS = ("slight_halophile", "moderate_halophile", "extreme_halophile")


def load_existing_ids() -> set[str]:
    rows = read_tsv(CURRENT_SHORTLIST)
    return {row["bacdive_id"] for row in rows}


def rank_candidate(row: Dict[str, str]) -> Tuple[int, int, str]:
    species = row["species"]
    genus = species.split()[0] if species else ""
    if base.EXTREME_GENUS_RE.search(genus):
        priority = 0
    elif base.MODERATE_GENUS_RE.search(genus):
        priority = 1
    else:
        priority = 2

    is_type = 0 if row.get("is_type_strain") == "yes" else 1
    return (priority, is_type, species.lower())


def fetch_and_classify(row: Dict[str, str]) -> Dict[str, object]:
    record = base.fetch_record(row["bacdive_id"])
    opt_entries = base.extract_optimum_halophily(record)
    chosen = None
    if opt_entries:
        opt_entries.sort(key=lambda e: len(base.parse_numbers(str(e.get("concentration", "")))))
        chosen = opt_entries[0]

    conc = base.parse_concentration_to_molar(chosen["concentration"]) if chosen else None
    label, include, reason = base.classify_halophile(conc)

    return {
        "record": record,
        "chosen_entry": chosen,
        "conc": conc,
        "label": label,
        "include": include,
        "reason": reason,
    }


def select_expansion_rows(
    rows: List[Dict[str, str]],
    per_class: int,
    require_type_strain: bool,
) -> List[Tuple[Dict[str, str], Dict[str, object]]]:
    existing_ids = load_existing_ids()
    selected: List[Tuple[Dict[str, str], Dict[str, object]]] = []
    counts: Counter[str] = Counter()

    candidate_rows = [row for row in rows if row["bacdive_id"] not in existing_ids]
    if require_type_strain:
        candidate_rows = [row for row in candidate_rows if row.get("is_type_strain") == "yes"]
    candidate_rows.sort(key=rank_candidate)

    for row in candidate_rows:
        if all(counts[label] >= per_class for label in TARGET_LABELS):
            break

        data = fetch_and_classify(row)
        label = str(data["label"] or "")
        if data["include"] != "yes":
            continue
        if label not in TARGET_LABELS:
            continue
        if counts[label] >= per_class:
            continue

        selected.append((row, data))
        counts[label] += 1

    return selected


def build_output_paths(batch_name: str) -> Dict[str, Path]:
    suffix = batch_name.strip()
    return {
        "screening": BASE / "curation" / "paper_screening_logs" / f"paper_screening_log_halophile_bacdive_{suffix}.tsv",
        "core": BASE / "curation" / "core_taxon_tables" / f"core_taxon_table_halophile_bacdive_{suffix}.tsv",
        "shortlist": BASE / "curation" / "source_tables" / f"bacdive_halophile_shortlist_{suffix}.tsv",
    }


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Generate a second-pass BacDive-backed halophile expansion batch."
    )
    parser.add_argument(
        "--per-class",
        type=int,
        default=5,
        help="Number of new clear-optimum rows to collect per halophile class.",
    )
    parser.add_argument(
        "--batch-name",
        default="expansion_round2",
        help="Suffix used in output filenames.",
    )
    parser.add_argument(
        "--allow-non-type-strains",
        action="store_true",
        help="Allow non-type strains into the expansion batch.",
    )
    args = parser.parse_args()

    # Workaround: the local Python had no SSL certificates for the BacDive API.
    ssl._create_default_https_context = ssl._create_unverified_context

    rows = base.read_candidates()
    selected = select_expansion_rows(
        rows=rows,
        per_class=args.per_class,
        require_type_strain=not args.allow_non_type_strains,
    )

    outputs = build_output_paths(args.batch_name)
    screening_rows = [base.screening_row(candidate, data) for candidate, data in selected]
    core_rows = [core_row_with_aliases(candidate, data) for candidate, data in selected]
    shortlist_rows = [base.shortlist_row(candidate, data) for candidate, data in selected]

    if screening_rows:
        write_tsv(outputs["screening"], screening_rows, list(screening_rows[0].keys()))
        write_tsv(outputs["core"], core_rows, list(core_rows[0].keys()))
        write_tsv(outputs["shortlist"], shortlist_rows, list(shortlist_rows[0].keys()))
    else:
        raise SystemExit("No expansion rows were selected.")

    class_counts = Counter(row["halophile_class"] for row in core_rows)
    print(
        json.dumps(
            {
                "selected_total": len(selected),
                "class_counts": dict(class_counts),
                "screening_out": str(outputs["screening"]),
                "core_out": str(outputs["core"]),
                "shortlist_out": str(outputs["shortlist"]),
            },
            indent=2,
        )
    )


def core_row_with_aliases(candidate_row: Dict[str, str], data: Dict[str, object]) -> Dict[str, str]:
    row = base.core_row(candidate_row, data)
    row["strain_aliases"] = candidate_row.get("strain_number_candidate", "")
    return row


if __name__ == "__main__":
    main()
