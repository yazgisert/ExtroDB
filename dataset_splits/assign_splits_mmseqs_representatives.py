"""
assign_splits_mmseqs_representatives.py

  1. Read the cluster representative sequences from mmseqs2_cluster_rep_seq.fasta
     (one protein per cluster at 30% identity — already deduplicated).
  2. Look up each representative's phenotype label via its proteome_id prefix.
  3. Randomly split representatives into train / val / test (70/10/20),
     stratified by phenotype label so every class is proportionally represented.
  4. Write output FASTAs for both splits:
       Split 1 (coarse, 5 classes):   rep_split1_train.faa / _val.faa / _test.faa
       Split 2 (subtypes, 12 classes): rep_split2_train.faa / _val.faa / _test.faa

Reproducibility
---------------
RANDOM_SEED = 42. Document this in the thesis methods section.
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

SPLITS_DIR   = Path(__file__).resolve().parent
META_TSV     = SPLITS_DIR / "split_input_metadata.tsv"
REP_FASTA    = SPLITS_DIR / "mmseqs2_cluster_rep_seq.fasta"

PARTITIONS   = ["train", "val", "test"]

# ── Helpers ──────────────────────────────────────────────────────────────────

def read_metadata(path: Path) -> dict[str, dict]:
    """Returns dict: proteome_id -> {split1_label, split2_label, protein_count}"""
    meta = {}
    with path.open(encoding="utf-8", newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            meta[row["proteome_id"]] = {
                "split1_label": row["phenotype_focus"],
                "split2_label": row["split2_label"],
            }
    return meta


def parse_representatives(fasta_path: Path, meta: dict) -> list[dict]:
    """
    Stream through the representative FASTA and collect one record per protein.
    Returns a list of dicts with keys: header, proteome_id, split1_label, split2_label.
    Proteins whose proteome_id is not in meta (shouldn't happen) are skipped.
    """
    records = []
    skipped = 0
    current_header = None

    with fasta_path.open(encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                current_header = line[1:]   # strip leading >
                proteome_id    = current_header.split("||")[0]
                info = meta.get(proteome_id)
                if info is None:
                    skipped += 1
                    current_header = None
                    continue
                records.append({
                    "header":       current_header,
                    "proteome_id":  proteome_id,
                    "split1_label": info["split1_label"],
                    "split2_label": info["split2_label"],
                    "seq_lines":    [],
                })
            elif current_header is not None and records:
                records[-1]["seq_lines"].append(line)

    if skipped:
        print(f"  WARNING: {skipped} representative(s) had unknown proteome_id — skipped.")
    return records


def stratified_split(
    records: list[dict],
    label_key: str,
    rng: random.Random,
) -> dict[str, str]:
    """
    Returns dict: record header -> partition name.
    Stratification is per label class.
    """
    by_label: dict[str, list[dict]] = defaultdict(list)
    for rec in records:
        by_label[rec[label_key]].append(rec)

    assignments: dict[str, str] = {}

    for label, group in sorted(by_label.items()):
        rng.shuffle(group)
        n = len(group)

        n_train = max(1, round(n * TRAIN_FRAC))
        n_val   = max(1, round(n * VAL_FRAC))
        n_test  = n - n_train - n_val

        if n_test < 1:
            n_val  = max(1, n_val - 1)
            n_test = n - n_train - n_val

        if n < 10:
            print(f"  NOTE: {label} has only {n} representatives → "
                  f"train={n_train}, val={n_val}, test={max(0,n_test)}")

        sizes = {"train": n_train, "val": n_val, "test": max(0, n_test)}
        idx = 0
        for partition in PARTITIONS:
            for rec in group[idx : idx + sizes[partition]]:
                assignments[rec["header"]] = partition
            idx += sizes[partition]

    return assignments


def write_fastas(
    records: list[dict],
    assignments: dict[str, str],
    label_key: str,
    split_name: str,
) -> dict[str, int]:
    """Write one FASTA per partition. Returns protein counts per partition."""
    handles = {
        p: (SPLITS_DIR / f"rep_{split_name}_{p}.faa").open("w", encoding="utf-8")
        for p in PARTITIONS
    }
    counts = {p: 0 for p in PARTITIONS}

    for rec in records:
        header    = rec["header"]
        partition = assignments.get(header)
        if partition is None:
            continue
        label = rec[label_key]
        fh    = handles[partition]
        fh.write(f">{header}||{label}\n")
        for seq_line in rec["seq_lines"]:
            fh.write(seq_line + "\n")
        counts[partition] += 1

    for fh in handles.values():
        fh.close()

    return counts


def write_summary(
    records: list[dict],
    assignments: dict[str, str],
    label_key: str,
    out_path: Path,
) -> None:
    counts: dict = defaultdict(lambda: defaultdict(int))
    for rec in records:
        label     = rec[label_key]
        partition = assignments.get(rec["header"], "unassigned")
        counts[label][partition] += 1

    with out_path.open("w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["label", "partition", "n_representatives", "pct"])
        for label in sorted(counts):
            total = sum(counts[label][p] for p in PARTITIONS)
            for partition in PARTITIONS:
                n   = counts[label][partition]
                pct = n / total * 100 if total > 0 else 0.0
                w.writerow([label, partition, n, f"{pct:.1f}%"])
    print(f"  Written: {out_path.name}")


# ── Main ─────────────────────────────────────────────────────────────────────

def main() -> None:
    print(f"Reading metadata: {META_TSV.name}")
    meta = read_metadata(META_TSV)
    print(f"  {len(meta)} organisms loaded")

    print(f"\nParsing representatives: {REP_FASTA.name}")
    records = parse_representatives(REP_FASTA, meta)
    print(f"  {len(records):,} representative proteins loaded")

    print(f"\n── Split 1: coarse phenotype (5 classes) ──")
    rng = random.Random(RANDOM_SEED)
    s1  = stratified_split(records, "split1_label", rng)
    counts1 = write_fastas(records, s1, "split1_label", "split1")
    for p, n in counts1.items():
        print(f"  rep_split1_{p}.faa: {n:,} proteins")
    write_summary(records, s1, "split1_label",
                  SPLITS_DIR / "rep_split1_summary.tsv")

    print(f"\n── Split 2: fine-grained subtypes ──")
    rng = random.Random(RANDOM_SEED)
    s2  = stratified_split(records, "split2_label", rng)
    counts2 = write_fastas(records, s2, "split2_label", "split2")
    for p, n in counts2.items():
        print(f"  rep_split2_{p}.faa: {n:,} proteins")
    write_summary(records, s2, "split2_label",
                  SPLITS_DIR / "rep_split2_summary.tsv")

    print(f"\nTotal representatives: {len(records):,}")
    print(f"Split 1 total written: {sum(counts1.values()):,}")
    print(f"Split 2 total written: {sum(counts2.values()):,}")
    print("\nDone.")


if __name__ == "__main__":
    main()
