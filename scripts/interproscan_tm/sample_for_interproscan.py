# Samples 500 proteins per class from the prokaryote splits and the algae proteomes for InterProScan.

import csv
import random
from collections import defaultdict
from pathlib import Path

SEED      = 42
N_PER_CLASS = 500
random.seed(SEED)

BASE_DIR  = Path(__file__).resolve().parents[2]
OUT_DIR   = BASE_DIR / "results" / "interproscan"
OUT_DIR.mkdir(parents=True, exist_ok=True)

CLASS_NAMES = ["acidophile", "alkaliphile", "halophile", "psychrophile", "thermophile"]

def parse_prokaryote_faas():
    faa_files = [
        BASE_DIR / "dataset_splits" / "rep_split1_train.faa",
        BASE_DIR / "dataset_splits" / "rep_split1_val.faa",
        BASE_DIR / "dataset_splits" / "rep_split1_test.faa",
    ]

    by_class = defaultdict(list)

    for faa in faa_files:
        current_id  = None
        current_cls = None
        current_seq = []

        with open(faa) as f:
            for line in f:
                line = line.strip()
                if line.startswith(">"):
                    if current_id and current_cls:
                        by_class[current_cls].append((current_id, "".join(current_seq)))
                    parts = line[1:].split("||")
                    current_id  = parts[1].split()[0].strip()
                    current_cls = parts[2].strip()
                    current_seq = []
                else:
                    current_seq.append(line)

        if current_id and current_cls:
            by_class[current_cls].append((current_id, "".join(current_seq)))

    return by_class


def parse_eukaryote_faas():
    label_tsv = BASE_DIR / "embeddings" / "eukaryote_algae_labels.tsv"
    index_tsv = BASE_DIR / "searchable_resource" / "extremophile_master_index.tsv"

    pid2class = {}
    with open(label_tsv) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            pid2class[row["protein_id"]] = row["broad_class"]

    fasta_paths = []
    with open(index_tsv) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            if row.get("domain_scope") != "eukaryote":
                continue
            fp = (row.get("local_fasta_path") or "").strip()
            if fp:
                fasta_paths.append(BASE_DIR / fp)

    by_class = defaultdict(list)

    for fasta_path in fasta_paths:
        if not fasta_path.exists():
            print(f"  WARNING: not found: {fasta_path}")
            continue

        current_id  = None
        current_seq = []

        with open(fasta_path) as f:
            for line in f:
                line = line.strip()
                if line.startswith(">"):
                    if current_id and current_id in pid2class:
                        cls = pid2class[current_id]
                        by_class[cls].append((current_id, "".join(current_seq)))
                    current_id  = line[1:].split()[0]
                    current_seq = []
                else:
                    current_seq.append(line)

        if current_id and current_id in pid2class:
            cls = pid2class[current_id]
            by_class[cls].append((current_id, "".join(current_seq)))

    return by_class


def sample_and_write(by_class, out_faa, label, n=N_PER_CLASS):
    rng = random.Random(SEED)
    sampled = {}

    for cls in CLASS_NAMES:
        seqs = by_class.get(cls, [])
        if not seqs:
            print(f"  [{label}] {cls}: no sequences — skipping")
            continue
        k = min(n, len(seqs))
        sampled[cls] = rng.sample(seqs, k)
        print(f"  [{label}] {cls}: {k:,} sequences sampled (from {len(seqs):,})")

    with open(out_faa, "w") as f:
        for cls, seqs in sampled.items():
            for pid, seq in seqs:
                f.write(f">{pid}||{cls}\n{seq}\n")

    total = sum(len(v) for v in sampled.values())
    print(f"  Written {total:,} sequences to {out_faa}\n")
    return sampled


if __name__ == "__main__":
    print("Parsing prokaryote FAA files ...")
    prok = parse_prokaryote_faas()
    for cls in CLASS_NAMES:
        print(f"  {cls}: {len(prok.get(cls, [])):,}")

    print("\nSampling prokaryotes ...")
    sample_and_write(prok, OUT_DIR / "prokaryote_sample.faa", "prokaryote")

    print("Parsing eukaryote FAA files ...")
    euk = parse_eukaryote_faas()
    for cls in CLASS_NAMES:
        print(f"  {cls}: {len(euk.get(cls, [])):,}")

    print("\nSampling eukaryotes ...")
    sample_and_write(euk, OUT_DIR / "eukaryote_sample.faa", "eukaryote")

    print("Done.")
