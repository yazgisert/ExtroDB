# Checks that the current biocentral server gives the same ESM2 650M embeddings as the stored ones (computed on
# the old server): re-embeds 20 random test proteins and compares them.

import random

import h5py
import numpy as np
from biocentral_api import BiocentralAPI, CommonEmbedder

N_PROTEINS = 20
STORED_H5  = "embeddings/split1_30pct_esm2.h5"
SPLIT_FAA  = "dataset_splits/rep_split1_test.faa"


def read_faa(filepath):
    seqs = {}
    current_id, current_seq = None, []
    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if current_id is not None:
                    seqs[current_id] = "".join(current_seq)
                current_id, current_seq = line[1:].split("||")[1].split()[0].strip(), []
            else:
                current_seq.append(line)
    if current_id is not None:
        seqs[current_id] = "".join(current_seq)
    return seqs


if __name__ == "__main__":
    seqs = read_faa(SPLIT_FAA)
    with h5py.File(STORED_H5, "r") as f:
        ids = random.Random(0).sample([i for i in seqs if i in f], N_PROTEINS)
        stored = {i: f[i][:] for i in ids}

    api = BiocentralAPI().wait_until_healthy()
    fresh = api.embed(
        embedder_name=CommonEmbedder.ESM2_650M,
        reduce=True,
        sequence_data={i: seqs[i] for i in ids},
        use_half_precision=False,
    ).run_with_progress().to_dict()

    max_diff = max(np.abs(stored[i] - fresh[i]).max() for i in ids)
    min_cos  = min(np.dot(stored[i], fresh[i]) / (np.linalg.norm(stored[i]) * np.linalg.norm(fresh[i])) for i in ids)
    print(f"\n{N_PROTEINS} proteins compared (old server vs. current server, ESM2 650M)")
    print(f"  largest absolute difference: {max_diff:.2e}")
    print(f"  lowest cosine similarity:    {min_cos:.6f}")
    print("  -> consistent" if max_diff < 1e-3 else "  -> different from the stored embeddings")
