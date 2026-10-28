# Clusters the embedded non-extremophile proteins at 30% identity (MMseqs2 easy-cluster, default coverage 80%);
# identical sequences are clustered once. The representatives are the input for make_clustered_binary_split.py.

import argparse
import subprocess
import sys
from collections import Counter, defaultdict
from pathlib import Path

import h5py

NONEXTREM_H5  = "embeddings/nonextremophile.h5"
NONEXTREM_FAA = "nonextremophile_data/nonextremophile_filtered.faa"
OUT_DIR       = "nonextremophile_data/clustering_30pct"

MIN_SEQ_ID = 0.3
NEEDED     = 99_964 + 14_289 + 28_576   # extremophile train + val + test in split1 30%


def read_faa(filepath):
    seqs = {}
    current_id, current_seq = None, []
    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if current_id is not None:
                    seqs[current_id] = "".join(current_seq)
                current_id, current_seq = line[1:].split()[0], []
            else:
                current_seq.append(line)
    if current_id is not None:
        seqs[current_id] = "".join(current_seq)
    return seqs


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--faa", default=NONEXTREM_FAA)
    parser.add_argument("--out-dir", default=OUT_DIR)
    args = parser.parse_args()
    out_dir    = Path(args.out_dir)
    unique_faa = out_dir / "nonextremophile_unique.faa"
    id_map_tsv = out_dir / "unique_to_all_ids.tsv"
    prefix     = out_dir / "nonext_cluster"
    clusters   = out_dir / "nonext_cluster_cluster.tsv"

    if clusters.exists():
        sys.exit(f"{clusters} already exists. Delete {out_dir} to cluster again.")
    out_dir.mkdir(parents=True, exist_ok=True)

    seqs = read_faa(args.faa)
    with h5py.File(NONEXTREM_H5, "r") as f:
        embedded = [protein_id for protein_id in f.keys() if protein_id in seqs]

    # Identical sequences are clustered once; every ID sharing the sequence follows its representative
    ids_by_seq = defaultdict(list)
    for protein_id in embedded:
        ids_by_seq[seqs[protein_id]].append(protein_id)

    with open(unique_faa, "w") as faa, open(id_map_tsv, "w") as id_map:
        id_map.write("unique_id\tall_ids\n")
        for seq, ids in ids_by_seq.items():
            faa.write(f">{ids[0]}\n{seq}\n")
            id_map.write(f"{ids[0]}\t{','.join(ids)}\n")
    print(f"{len(embedded):,} embedded non-extremophile proteins, {len(ids_by_seq):,} unique sequences")

    cmd = ["mmseqs", "easy-cluster", str(unique_faa), str(prefix), str(out_dir / "mmseqs_tmp"),
           "--min-seq-id", str(MIN_SEQ_ID), "--threads", "8"]
    print("\nRunning MMseqs2 clustering at 30% identity (default coverage 80%) ...")
    print("  " + " ".join(cmd))
    subprocess.run(cmd, check=True)

    cluster_of = {}
    with open(clusters) as f:
        for line in f:
            rep, member = line.rstrip("\n").split("\t")
            cluster_of[member] = rep
    sizes = Counter(cluster_of.values())
    n_reps = len(sizes)

    lines = [
        f"Unique sequences clustered:  {len(cluster_of):,}",
        f"Clusters (representatives):  {n_reps:,}",
        f"Singletons:                  {sum(s == 1 for s in sizes.values()):,}",
        f"Largest cluster:             {max(sizes.values()):,}",
        f"Needed for a 50/50 split:    {NEEDED:,}  (99,964 train + 14,289 val + 28,576 test)",
        f"Representatives vs. needed:  {n_reps / NEEDED:.1%}",
    ]
    summary = "\n".join(lines)
    print("\n" + summary)
    (out_dir / "cluster_summary.txt").write_text(summary + "\n")
    print(f"\nSaved to {out_dir}")
