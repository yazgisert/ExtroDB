#!/usr/bin/env python3

from __future__ import annotations

import sys
from pathlib import Path
from typing import Dict, List

from common import read_tsv, write_tsv


COMMON_FIELDS = {
    "taxon_name": "taxon_name",
    "ncbi_taxid": "ncbi_taxid",
    "strain": "strain",
    "strain_aliases": "strain_aliases",
    "evidence_source": "evidence_source",
}
MODE_FIELDS = {
    "psychrophile": ["psychrophile", "ogt_celsius"],
    "halophile": ["halophile", "halophile_class", "nacl_optimum_M"],
    "thermophile": ["thermophile", "hyperthermophile", "ogt_celsius"],
    "acidophile": ["acidophile", "acidophile_class", "ph_optimum"],
    "acidotolerant": ["acidophile", "acidophile_class", "ph_optimum"],
    "alkaliphile": ["alkaliphile", "alkaliphile_class", "ph_optimum"],
}
TRACKER_DEFAULTS = {
    "sequence_source_db": "",
    "sequence_source_accession": "",
    "sequence_level": "",
    "retrieval_status": "pending",
    "retrieval_notes": "",
}


def infer_mode(rows: List[Dict[str, str]], core_path: Path) -> str:
    # Prefer an explicit mode hint from the input filename because the shared
    # core template can now carry multiple phenotype column groups at once.
    stem = core_path.stem.lower()
    for mode in ("acidotolerant", "acidophile", "alkaliphile", "halophile", "psychrophile", "thermophile"):
        if mode in stem:
            return mode

    # Fall back to phenotype columns for older files that predate the shared
    # all-phenotypes template.
    first = rows[0]
    if "psychrophile" in first and "ogt_celsius" in first:
        return "psychrophile"
    if "acidophile" in first and "acidophile_class" in first and "ph_optimum" in first:
        return "acidophile"
    if "alkaliphile" in first and "alkaliphile_class" in first and "ph_optimum" in first:
        return "alkaliphile"
    if "halophile" in first and "halophile_class" in first and "nacl_optimum_M" in first:
        return "halophile"
    if "thermophile" in first and "hyperthermophile" in first and "ogt_celsius" in first:
        return "thermophile"
    raise ValueError("Could not infer tracker mode from input columns.")


def build_tracker_rows(rows: List[Dict[str, str]], mode: str) -> List[Dict[str, str]]:
    # A tracker starts as a copy of the curated core table plus blank retrieval
    # fields that will be filled in later during sequence collection.
    return [
        {
            **{field: row.get(source, "") for field, source in COMMON_FIELDS.items()},
            **{field: row.get(field, "") for field in MODE_FIELDS[mode]},
            **TRACKER_DEFAULTS,
        }
        for row in rows
    ]


def main() -> int:
    if len(sys.argv) != 3:
        print("Usage: initialize_sequence_tracker.py CORE_TSV TRACKER_TSV", file=sys.stderr)
        return 1

    core_path = Path(sys.argv[1])
    tracker_path = Path(sys.argv[2])

    rows = read_tsv(core_path)
    if not rows:
        raise SystemExit("Input core TSV has no rows.")

    mode = infer_mode(rows, core_path)
    tracker_rows = build_tracker_rows(rows, mode)
    write_tsv(tracker_path, tracker_rows, list(tracker_rows[0].keys()))

    print(f"mode={mode}")
    print(f"rows={len(tracker_rows)}")
    print(f"tracker={tracker_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
