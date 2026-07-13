#!/usr/bin/env python3
"""
Starter script for joining AFESM-style taxonomy and biome tables.

This script is intentionally conservative:
- it assumes you have already downloaded the processed AFESM outputs
- it joins on representative ID
- it adds simple first-pass flags for biome relevance and taxonomy specificity
- it does NOT decide final extremophile labels automatically

Expected AFESM-style inputs, based on the afesm-analysis repository:

1. Taxonomy table
   Example output name:
   afesm30_repseq_foldseek_clu_nonsingleton_allMem_lca_blacklist-
   repId_taxId_taxRank_taxName_groupName.tsv

   Expected columns in order:
   rep_id, tax_id, tax_rank, tax_name, group_name

2. Biome table
   Example output name:
   afesm30-repId_nBiome_nLcaLen_biomeLca.tsv

   Expected columns in order:
   rep_id, n_biomes, biome_lca_depth, biome_lca

Output:
   A joined candidate-discovery TSV suitable for later curation.
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path


def load_tsv(path: Path, columns: list[str]) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        for line_number, row in enumerate(reader, start=1):
            if not row:
                continue
            if len(row) < len(columns):
                raise ValueError(
                    f"{path} line {line_number}: expected at least {len(columns)} "
                    f"tab-separated fields, found {len(row)}"
                )
            rows.append({columns[i]: row[i].strip() for i in range(len(columns))})
    return rows


def normalize_taxonomy_specificity(tax_rank: str) -> str:
    rank = tax_rank.strip().lower()
    if rank in {"species", "subspecies", "strain", "genus"}:
        return "high"
    if rank in {"family", "order"}:
        return "medium"
    if rank in {"class", "phylum", "kingdom", "superkingdom", "root", "domain", ""}:
        return "low"
    return "unknown"


def normalize_biome_category(biome_lca: str) -> str:
    biome = biome_lca.strip().lower()
    if "salt" in biome:
        return "salt"
    if "marine" in biome:
        return "marine"
    if "soil" in biome:
        return "soil"
    if "wastewater" in biome:
        return "wastewater"
    if "human" in biome or "host" in biome:
        return "host"
    if "lab" in biome or "engineered" in biome:
        return "engineered"
    if "air" in biome:
        return "air"
    if "food" in biome:
        return "food"
    if "plant" in biome:
        return "plant"
    return "unknown"


def infer_extreme_biome_flag(category: str) -> str:
    return "yes" if category == "salt" else "no"


def infer_candidate_type(category: str) -> str:
    if category == "salt":
        return "halophile_candidate"
    if category == "marine":
        return "marine_candidate"
    return "unclassified"


def build_joined_rows(
    taxonomy_rows: list[dict[str, str]],
    biome_rows: list[dict[str, str]],
    source_label: str,
) -> list[dict[str, str]]:
    biome_index = {row["rep_id"]: row for row in biome_rows}
    joined: list[dict[str, str]] = []

    for tax in taxonomy_rows:
        rep_id = tax["rep_id"]
        biome = biome_index.get(rep_id, {})
        biome_lca = biome.get("biome_lca", "")
        biome_category = normalize_biome_category(biome_lca)

        joined.append(
            {
                "afesm_rep_id": rep_id,
                "predicted_tax_id": tax.get("tax_id", ""),
                "predicted_tax_name": tax.get("tax_name", ""),
                "predicted_tax_rank": tax.get("tax_rank", ""),
                "taxonomy_group_name": tax.get("group_name", ""),
                "taxonomy_specificity_tier": normalize_taxonomy_specificity(
                    tax.get("tax_rank", "")
                ),
                "biome_lca": biome_lca,
                "n_biomes": biome.get("n_biomes", ""),
                "biome_lca_depth": biome.get("biome_lca_depth", ""),
                "extreme_biome_category": biome_category,
                "extreme_biome_flag": infer_extreme_biome_flag(biome_category),
                "extremophile_candidate_type": infer_candidate_type(biome_category),
                "candidate_source": source_label,
                "evidence_type": "afesm_taxonomy_biome_only",
                "include_in_core_dataset": "review",
                "review_notes": "",
            }
        )

    return joined


def filter_rows(
    rows: list[dict[str, str]],
    biome_keyword: str | None,
    min_n_biomes: int | None,
    require_tax_specificity: bool,
) -> list[dict[str, str]]:
    filtered = rows

    if biome_keyword:
        keyword = biome_keyword.lower()
        filtered = [r for r in filtered if keyword in r["biome_lca"].lower()]

    if min_n_biomes is not None:
        def has_min_biomes(row: dict[str, str]) -> bool:
            try:
                return int(row["n_biomes"]) >= min_n_biomes
            except ValueError:
                return False

        filtered = [r for r in filtered if has_min_biomes(r)]

    if require_tax_specificity:
        filtered = [
            r
            for r in filtered
            if r["taxonomy_specificity_tier"] in {"high", "medium"}
        ]

    return filtered


def write_tsv(path: Path, rows: list[dict[str, str]]) -> None:
    if not rows:
        raise ValueError("No rows to write after filtering.")

    fieldnames = list(rows[0].keys())
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Join AFESM-style taxonomy and biome tables into a candidate TSV."
    )
    parser.add_argument(
        "--taxonomy-tsv",
        required=True,
        type=Path,
        help="Path to the AFESM-style taxonomy TSV.",
    )
    parser.add_argument(
        "--biome-tsv",
        required=True,
        type=Path,
        help="Path to the AFESM-style biome TSV.",
    )
    parser.add_argument(
        "--output-tsv",
        required=True,
        type=Path,
        help="Path for the joined output TSV.",
    )
    parser.add_argument(
        "--source-label",
        default="Yeo_2026_AFESM",
        help="Provenance label written into candidate_source.",
    )
    parser.add_argument(
        "--biome-keyword",
        default=None,
        help="Optional substring filter for biome_lca, e.g. Salt or Marine.",
    )
    parser.add_argument(
        "--min-n-biomes",
        type=int,
        default=None,
        help="Optional minimum n_biomes threshold.",
    )
    parser.add_argument(
        "--require-tax-specificity",
        action="store_true",
        help="Keep only rows with medium/high taxonomy specificity.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()

    taxonomy_columns = ["rep_id", "tax_id", "tax_rank", "tax_name", "group_name"]
    biome_columns = ["rep_id", "n_biomes", "biome_lca_depth", "biome_lca"]

    taxonomy_rows = load_tsv(args.taxonomy_tsv, taxonomy_columns)
    biome_rows = load_tsv(args.biome_tsv, biome_columns)

    joined_rows = build_joined_rows(
        taxonomy_rows=taxonomy_rows,
        biome_rows=biome_rows,
        source_label=args.source_label,
    )

    filtered_rows = filter_rows(
        rows=joined_rows,
        biome_keyword=args.biome_keyword,
        min_n_biomes=args.min_n_biomes,
        require_tax_specificity=args.require_tax_specificity,
    )

    write_tsv(args.output_tsv, filtered_rows)

    print(f"Wrote {len(filtered_rows)} joined candidate rows to {args.output_tsv}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
