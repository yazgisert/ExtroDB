# Splits the MMseqs2 cluster representatives (30% identity by default; --rep-fasta / --tag for the 20% run) into
# train/val/test (70/10/20, stratified per class, seed 42), downsamples every class to the smallest one within
# each partition, and writes the FASTAs and summaries for split1 (5 classes) and split2 (12 classes).

from __future__ import annotations

import argparse
import csv
import random
from collections import defaultdict
from pathlib import Path


RANDOM_SEED = 42
TRAIN_FRAC  = 0.70
VAL_FRAC    = 0.10
TEST_FRAC   = 0.20

SPLITS_DIR   = Path(__file__).resolve().parents[2] / "dataset_splits"
META_TSV     = SPLITS_DIR / "split_input_metadata.tsv"
REP_FASTA    = SPLITS_DIR / "mmseqs2_cluster_rep_seq.fasta"

PARTITIONS   = ["train", "val", "test"]


def read_metadata(path: Path) -> dict[str, dict]:
    meta = {}
    with path.open(encoding="utf-8", newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            meta[row["proteome_id"]] = {
                "split1_label": row["phenotype_focus"],
                "split2_label": row["split2_label"],
            }
    return meta


def parse_representatives(fasta_path: Path, meta: dict) -> list[dict]:
    records = []
    skipped = 0
    current_header = None

    with fasta_path.open(encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                current_header = line[1:]
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


def balance_assignments(
    records: list[dict],
    assignments: dict[str, str],
    label_key: str,
    rng: random.Random,
) -> dict[str, str]:
    by_label_partition: dict[str, dict[str, list[str]]] = defaultdict(
        lambda: defaultdict(list)
    )
    for rec in records:
        header    = rec["header"]
        partition = assignments.get(header)
        if partition is None:
            continue
        by_label_partition[rec[label_key]][partition].append(header)

    class_totals = {
        label: sum(len(v) for v in partitions.values())
        for label, partitions in by_label_partition.items()
    }
    min_label = min(class_totals, key=class_totals.get)
    min_size  = class_totals[min_label]
    print(f"  Smallest class: '{min_label}' ({min_size:,} representatives) "
          f"— downsampling all larger classes to match.")

    balanced: dict[str, str] = {}
    for label, partitions in sorted(by_label_partition.items()):
        total = class_totals[label]
        for partition, headers in partitions.items():
            if total <= min_size:
                # Already at or below the cap — keep everything
                keep = headers
            else:
                # Keep the proportional share of min_size from this partition
                n_keep = round(min_size * len(headers) / total)
                keep   = rng.sample(headers, min(n_keep, len(headers)))
            for h in keep:
                balanced[h] = partition

        if total > min_size:
            kept = sum(
                round(min_size * len(headers) / total)
                for headers in partitions.values()
            )
            print(f"  {label}: {total:,} → ~{kept:,} representatives")

    return balanced


def write_fastas(
    records: list[dict],
    assignments: dict[str, str],
    label_key: str,
    split_name: str,
    tag: str = "",
) -> dict[str, int]:
    handles = {
        p: (SPLITS_DIR / f"rep_{split_name}{tag}_{p}.faa").open("w", encoding="utf-8")
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
    tag: str = "",
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


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Assign MMseqs2 cluster representatives to train/val/test splits."
    )
    p.add_argument(
        "--rep-fasta",
        default=None,
        help=(
            "Path to the representative FASTA produced by mmseqs2 easy-cluster. "
            f"Default: {REP_FASTA.name}"
        ),
    )
    p.add_argument(
        "--tag",
        default="",
        help=(
            "Suffix added to output file names before the partition, e.g. "
            "'_20pct' → rep_split1_20pct_train.faa. Default: '' (no suffix)."
        ),
    )
    return p.parse_args()


def main() -> None:
    args     = parse_args()
    rep_fasta = Path(args.rep_fasta) if args.rep_fasta else REP_FASTA
    tag       = args.tag

    print(f"Reading metadata: {META_TSV.name}")
    meta = read_metadata(META_TSV)
    print(f"  {len(meta)} organisms loaded")

    print(f"\nParsing representatives: {rep_fasta.name}")
    records = parse_representatives(rep_fasta, meta)
    print(f"  {len(records):,} representative proteins loaded")

    print("\nSplit 1: coarse phenotype (5 classes)")
    rng = random.Random(RANDOM_SEED)
    s1  = stratified_split(records, "split1_label", rng)
    print(f"  Balancing classes:")
    s1  = balance_assignments(records, s1, "split1_label", rng)
    counts1 = write_fastas(records, s1, "split1_label", "split1", tag=tag)
    for p, n in counts1.items():
        print(f"  rep_split1{tag}_{p}.faa: {n:,} proteins")
    write_summary(records, s1, "split1_label",
                  SPLITS_DIR / f"rep_split1{tag}_summary.tsv")

    print("\nSplit 2: fine-grained subtypes (12 classes)")
    rng = random.Random(RANDOM_SEED)
    s2  = stratified_split(records, "split2_label", rng)
    print(f"  Balancing classes:")
    s2  = balance_assignments(records, s2, "split2_label", rng)
    counts2 = write_fastas(records, s2, "split2_label", "split2", tag=tag)
    for p, n in counts2.items():
        print(f"  rep_split2{tag}_{p}.faa: {n:,} proteins")
    write_summary(records, s2, "split2_label",
                  SPLITS_DIR / f"rep_split2{tag}_summary.tsv")

    print(f"\nTotal representatives: {len(records):,}")
    print(f"Split 1 total written: {sum(counts1.values()):,}")
    print(f"Split 2 total written: {sum(counts2.values()):,}")
    print("\nDone.")


if __name__ == "__main__":
    main()
