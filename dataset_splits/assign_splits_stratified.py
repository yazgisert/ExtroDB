"""
assign_splits_stratified.py

Assigns the 192 single-label prokaryotic organisms to train / val / test
using random stratified sampling — the 70/10/20 ratio is applied
independently within each phenotype class so every class is proportionally
represented in all three partitions.

Split 1  —  coarse label    (phenotype_focus: 5 classes)
Split 2  —  fine-grained    (split2_label: 12 subtypes)

Outputs (written alongside this script in dataset_splits/)
----------------------------------------------------------
  split1_stratified_assignments.tsv
  split2_stratified_assignments.tsv
  split1_stratified_summary.tsv
  split2_stratified_summary.tsv

Reproducibility
---------------
Set RANDOM_SEED below. The default (42) is arbitrary; document whatever
value you use in your thesis methods section.
"""

from __future__ import annotations

import csv
import random
from collections import defaultdict
from pathlib import Path

# ── Config ───────────────────────────────────────────────────────────────────

RANDOM_SEED = 42
TRAIN_FRAC  = 0.70
VAL_FRAC    = 0.10
TEST_FRAC   = 0.20

SPLITS_DIR = Path(__file__).resolve().parent
META_TSV   = SPLITS_DIR / "split_input_metadata.tsv"

OUT_S1     = SPLITS_DIR / "split1_stratified_assignments.tsv"
OUT_S2     = SPLITS_DIR / "split2_stratified_assignments.tsv"
OUT_S1_SUM = SPLITS_DIR / "split1_stratified_summary.tsv"
OUT_S2_SUM = SPLITS_DIR / "split2_stratified_summary.tsv"

PARTITIONS = ["train", "val", "test"]

# ── Helpers ──────────────────────────────────────────────────────────────────

def read_metadata(path: Path) -> list[dict]:
    with path.open(encoding="utf-8", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def stratified_split(
    organisms: list[dict],
    label_key: str,
    rng: random.Random,
) -> dict[str, str]:
    """
    Returns dict: proteome_id -> partition

    For each label class:
      - shuffle organisms
      - assign first 70% to train, next 10% to val, remainder to test
      - for very small classes (< 5 organisms) every partition gets at
        least 1 organism where possible; a warning is printed
    """
    # Group by label
    by_label: dict[str, list[dict]] = defaultdict(list)
    for org in organisms:
        by_label[org[label_key]].append(org)

    assignments: dict[str, str] = {}

    for label, group in sorted(by_label.items()):
        rng.shuffle(group)
        n = len(group)

        if n < 3:
            # Too small to split meaningfully — all go to train
            print(f"  WARNING: {label} has only {n} organism(s) — all assigned to train")
            for org in group:
                assignments[org["proteome_id"]] = "train"
            continue

        n_train = max(1, round(n * TRAIN_FRAC))
        n_val   = max(1, round(n * VAL_FRAC))
        n_test  = n - n_train - n_val

        if n_test < 1:
            # Reduce val by 1 to free up a test slot
            n_val  = max(1, n_val - 1)
            n_test = n - n_train - n_val

        if n_test < 1:
            print(f"  WARNING: {label} has only {n} organism(s) — "
                  f"cannot fill all three partitions; test will be empty")
            n_test = 0

        sizes = {"train": n_train, "val": n_val, "test": n_test}

        if n < 10:
            print(f"  NOTE: {label} has only {n} organism(s) → "
                  f"train={n_train}, val={n_val}, test={n_test}")

        idx = 0
        for partition in PARTITIONS:
            for org in group[idx : idx + sizes[partition]]:
                assignments[org["proteome_id"]] = partition
            idx += sizes[partition]

    return assignments


def write_assignments(
    assignments: dict[str, str],
    meta_rows: list[dict],
    label_key: str,
    out_path: Path,
) -> None:
    meta = {r["proteome_id"]: r for r in meta_rows}
    with out_path.open("w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["proteome_id", "label", "partition", "protein_count"])
        for pid in sorted(assignments):
            w.writerow([
                pid,
                meta[pid][label_key],
                assignments[pid],
                meta[pid]["protein_count"],
            ])
    print(f"  Written: {out_path.name}")


def write_summary(
    assignments: dict[str, str],
    meta_rows: list[dict],
    label_key: str,
    out_path: Path,
) -> None:
    meta = {r["proteome_id"]: r for r in meta_rows}
    counts: dict = defaultdict(lambda: defaultdict(lambda: {"orgs": 0, "proteins": 0}))

    for pid, partition in assignments.items():
        label = meta[pid][label_key]
        counts[label][partition]["orgs"]     += 1
        counts[label][partition]["proteins"] += int(meta[pid]["protein_count"] or 0)

    with out_path.open("w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["label", "partition", "n_organisms", "n_proteins", "pct_proteins"])
        for label in sorted(counts):
            total = sum(counts[label][p]["proteins"] for p in PARTITIONS)
            for partition in PARTITIONS:
                d = counts[label][partition]
                pct = d["proteins"] / total * 100 if total > 0 else 0.0
                w.writerow([label, partition, d["orgs"], d["proteins"], f"{pct:.1f}%"])

    print(f"  Written: {out_path.name}")


def print_summary(path: Path) -> None:
    with path.open(encoding="utf-8") as fh:
        print(fh.read())


# ── Main ─────────────────────────────────────────────────────────────────────

def main() -> None:
    rng = random.Random(RANDOM_SEED)

    print(f"Reading metadata: {META_TSV.name}")
    meta_rows = read_metadata(META_TSV)
    print(f"  {len(meta_rows)} organisms loaded")
    print(f"  Random seed: {RANDOM_SEED}\n")

    print("── Split 1: coarse phenotype (5 classes) ──")
    s1 = stratified_split(meta_rows, "phenotype_focus", rng)
    write_assignments(s1, meta_rows, "phenotype_focus", OUT_S1)
    write_summary(s1, meta_rows, "phenotype_focus", OUT_S1_SUM)

    # Reset RNG so Split 2 is independently reproducible from the same seed
    rng = random.Random(RANDOM_SEED)

    print("\n── Split 2: fine-grained subtypes ──")
    s2 = stratified_split(meta_rows, "split2_label", rng)
    write_assignments(s2, meta_rows, "split2_label", OUT_S2)
    write_summary(s2, meta_rows, "split2_label", OUT_S2_SUM)

    print("\n── Split 1 summary ──")
    print_summary(OUT_S1_SUM)

    print("── Split 2 summary ──")
    print_summary(OUT_S2_SUM)

    print("Done.")
    print(f"\nAll outputs written to: {SPLITS_DIR}")


if __name__ == "__main__":
    main()
