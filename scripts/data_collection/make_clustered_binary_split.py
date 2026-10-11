# Builds the clustered binary and 6-class splits: the 13,031 non-extremophile cluster representatives (30%) are
# split 70/10/20, and each split gets the same number of extremophiles, equal per class, from the existing
# split1_30pct splits. 6-class: every class as large as the non-extremophile class in train/val; the test set is
# the binary one.

import sys
from collections import Counter, defaultdict
from pathlib import Path

import h5py
import numpy as np

SEED = 42

CLUSTERS  = "nonextremophile_data/clustering_30pct/nonext_cluster_cluster.tsv"
EXTREM_H5 = "embeddings/split1_30pct.h5"
EXTREM_FAA = {
    "train": "dataset_splits/rep_split1_train.faa",
    "val":   "dataset_splits/rep_split1_val.faa",
    "test":  "dataset_splits/rep_split1_test.faa",
}
CLASSES = ["acidophile", "alkaliphile", "halophile", "psychrophile", "thermophile"]
SPLITS  = ["train", "val", "test"]

OUT_DIR      = Path("dataset_splits/clustered_binary")
BINARY_TSV   = OUT_DIR / "binary_split.tsv"
SIXCLASS_TSV = OUT_DIR / "sixclass_split.tsv"


def extremophile_ids():
    # Extremophile representatives per split and class, restricted to proteins that have ProtT5 embeddings
    with h5py.File(EXTREM_H5, "r") as f:
        embedded = set(f.keys())
    ids = {split: defaultdict(list) for split in SPLITS}
    for split, path in EXTREM_FAA.items():
        with open(path) as f:
            for line in f:
                if line.startswith(">"):
                    parts = line[1:].split("||")
                    protein_id, label = parts[1].split()[0].strip(), parts[2].strip()
                    if protein_id in embedded:
                        ids[split][label].append(protein_id)
    for split in SPLITS:
        for label in CLASSES:
            ids[split][label].sort()
    return ids


def sample_per_class(pool, total, rng):
    # Same number from every extremophile class; the remainder goes to the first classes
    per_class = [total // len(CLASSES) + (i < total % len(CLASSES)) for i in range(len(CLASSES))]
    chosen = {}
    for label, n in zip(CLASSES, per_class):
        chosen[label] = list(rng.choice(pool[label], n, replace=False))
    return chosen


def write(path, rows):
    with open(path, "w") as f:
        f.write("protein_id\tsplit\tlabel\n")
        for protein_id, split, label in rows:
            f.write(f"{protein_id}\t{split}\t{label}\n")


if __name__ == "__main__":
    if BINARY_TSV.exists():
        sys.exit(f"{BINARY_TSV} already exists. Delete {OUT_DIR} to build the split again.")
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)

    with open(CLUSTERS) as f:
        reps = sorted({line.split("\t")[0] for line in f})
    rng.shuffle(reps)
    n_train, n_val = int(len(reps) * 0.7), int(len(reps) * 0.1)
    nonext = {"train": reps[:n_train], "val": reps[n_train:n_train + n_val], "test": reps[n_train + n_val:]}

    ext = extremophile_ids()

    binary_rows, binary_ext = [], {}
    for split in SPLITS:
        binary_ext[split] = sample_per_class(ext[split], len(nonext[split]), rng)
        binary_rows += [(i, split, "non_extremophile") for i in nonext[split]]
        binary_rows += [(i, split, label) for label, chosen in binary_ext[split].items() for i in chosen]
    write(BINARY_TSV, binary_rows)

    # 6-class: every extremophile class as large as the non-extremophile class; same test set as the binary model
    sixclass_rows = [row for row in binary_rows if row[1] == "test"]
    for split in ("train", "val"):
        sixclass_rows += [(i, split, "non_extremophile") for i in nonext[split]]
        for label in CLASSES:
            chosen = rng.choice(ext[split][label], len(nonext[split]), replace=False)
            sixclass_rows += [(i, split, label) for i in chosen]
    write(SIXCLASS_TSV, sixclass_rows)

    for name, rows in (("Binary", binary_rows), ("6-class", sixclass_rows)):
        counts = Counter((split, label) for _, split, label in rows)
        print(f"\n{name} split")
        print(f"  {'':18}" + "".join(f"{s:>9}" for s in SPLITS))
        for label in CLASSES + ["non_extremophile"]:
            print(f"  {label:<18}" + "".join(f"{counts[(s, label)]:>9,}" for s in SPLITS))
        print(f"  {'total':<18}" + "".join(f"{sum(v for (s, _), v in counts.items() if s == sp):>9,}" for sp in SPLITS))

    print(f"\nSaved to {BINARY_TSV}\n         {SIXCLASS_TSV}")
