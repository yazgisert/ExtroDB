# MMseqs2 search of the binary test proteins against train + val (--clustered for the clustered split). The hits feed
# the leakage check and the binary homology baseline.

import argparse
import subprocess
import sys
from pathlib import Path

import h5py
import numpy as np

SEED = 42

EXTREM_H5     = "embeddings/split1_30pct.h5"
NONEXTREM_H5  = "embeddings/nonextremophile.h5"
NONEXTREM_FAA = "nonextremophile_data/nonextremophile_filtered.faa"
EXTREM_FAA = {
    "train": "dataset_splits/rep_split1_train.faa",
    "val":   "dataset_splits/rep_split1_val.faa",
    "test":  "dataset_splits/rep_split1_test.faa",
}

CLUSTERED_SPLIT = "dataset_splits/clustered_binary/binary_split.tsv"

SENSITIVITY = 7.5
MAX_EVALUE  = 0.1


def read_faa(filepath, id_from_header):
    seqs, labels = {}, {}
    current_id, current_seq = None, []
    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if current_id is not None:
                    seqs[current_id] = "".join(current_seq)
                current_id, label = id_from_header(line[1:])
                labels[current_id] = label
                current_seq = []
            else:
                current_seq.append(line)
    if current_id is not None:
        seqs[current_id] = "".join(current_seq)
    return seqs, labels


def extremophile_header(header):
    parts = header.split("||")
    return parts[1].split()[0].strip(), parts[2].strip()


def nonextremophile_header(header):
    return header.split()[0], "non_extremophile"


def binary_split():
    # Original binary split (random non-extremophile split), with duplicated test sequences removed
    ext_seqs, ext_labels, ext_split = {}, {}, {}
    for split, path in EXTREM_FAA.items():
        seqs, labels = read_faa(path, extremophile_header)
        ext_seqs.update(seqs)
        ext_labels.update(labels)
        ext_split.update({i: split for i in seqs})

    members = {"train": [], "val": [], "test": []}
    with h5py.File(EXTREM_H5, "r") as f:
        for seq_id in f:
            if seq_id in ext_split:
                members[ext_split[seq_id]].append(seq_id)

    with h5py.File(NONEXTREM_H5, "r") as f:
        non_ids = list(f.keys())
    np.random.default_rng(SEED).shuffle(non_ids)
    n_train, n_val, n_test = (len(members[s]) for s in ("train", "val", "test"))
    non = {
        "train": non_ids[:n_train],
        "val":   non_ids[n_train:n_train + n_val],
        "test":  non_ids[n_train + n_val:n_train + n_val + n_test],
    }

    non_seqs, non_labels = read_faa(NONEXTREM_FAA, nonextremophile_header)
    seen = {non_seqs[i] for i in non["train"] + non["val"]}
    non["test"] = [i for i in non["test"] if non_seqs[i] not in seen]

    seqs   = {**ext_seqs, **non_seqs}
    labels = {**ext_labels, **non_labels}
    query     = members["test"] + non["test"]
    reference = members["train"] + members["val"] + non["train"] + non["val"]
    return query, reference, seqs, labels


def clustered_split():
    ext_seqs, ext_labels = {}, {}
    for path in EXTREM_FAA.values():
        seqs, labels = read_faa(path, extremophile_header)
        ext_seqs.update(seqs)
        ext_labels.update(labels)
    non_seqs, non_labels = read_faa(NONEXTREM_FAA, nonextremophile_header)
    seqs, labels = {**ext_seqs, **non_seqs}, {**ext_labels, **non_labels}

    query, reference = [], []
    with open(CLUSTERED_SPLIT) as f:
        next(f)
        for line in f:
            protein_id, split, _ = line.rstrip("\n").split("\t")
            (query if split == "test" else reference).append(protein_id)
    return query, reference, seqs, labels


def write_faa(path, ids, seqs):
    with open(path, "w") as f:
        for i in ids:
            f.write(f">{i}\n{seqs[i]}\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--clustered", action="store_true", help="use the clustered binary split instead of the original one")
    clustered = parser.parse_args().clustered

    OUT_DIR   = Path("results/test_vs_train_search_clustered" if clustered else "results/test_vs_train_search")
    QUERY_FAA = OUT_DIR / "test_query.faa"
    REF_FAA   = OUT_DIR / "train_reference.faa"
    LABELS    = OUT_DIR / "labels.tsv"
    HITS      = OUT_DIR / "hits.tsv"
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    print("Building query (test) and reference (train + val) ...")
    query, reference, seqs, labels = clustered_split() if clustered else binary_split()
    write_faa(QUERY_FAA, query, seqs)
    write_faa(REF_FAA, reference, seqs)
    with open(LABELS, "w") as f:
        f.write("protein_id\tset\tlabel\n")
        for i in query:
            f.write(f"{i}\ttest\t{labels[i]}\n")
        for i in reference:
            f.write(f"{i}\ttrain\t{labels[i]}\n")
    print(f"  query:     {len(query):,} test proteins")
    print(f"  reference: {len(reference):,} training proteins")

    if HITS.exists():
        sys.exit(f"\n{HITS} already exists. Delete it to search again.")

    cmd = [
        "mmseqs", "easy-search", str(QUERY_FAA), str(REF_FAA), str(HITS), str(OUT_DIR / "mmseqs_tmp"),
        "-s", str(SENSITIVITY),
        "-e", str(MAX_EVALUE),
        "--format-output", "query,target,pident,evalue,bits,qcov,tcov",
        "--threads", "8",
    ]
    print("\nRunning MMseqs2 (no identity cutoff, e-value <= 0.1) ...")
    print("  " + " ".join(cmd))
    subprocess.run(cmd, check=True)
    print(f"\nHits saved to {HITS}")
