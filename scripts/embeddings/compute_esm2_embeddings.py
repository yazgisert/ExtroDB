# ESM2 650M embeddings for the 5-class splits via the biocentral API. Batches of 1000, falling back to 100 and
# then 10 if a batch fails; resumes from the existing .h5 files.
# Run: caffeinate -i python scripts/embeddings/compute_esm2_embeddings.py

import h5py
import numpy as np
import os
import time
from biocentral_api import BiocentralAPI, CommonEmbedder, batched

ESM2_650M  = CommonEmbedder('facebook/esm2_t33_650M_UR50D')
BATCH_SIZE = 1000
OUT_DIR    = "embeddings"

SPLITS = {
    "split1_30pct": [
        "dataset_splits/rep_split1_train.faa",
        "dataset_splits/rep_split1_val.faa",
        "dataset_splits/rep_split1_test.faa",
    ],
    "split2_30pct": [
        "dataset_splits/rep_split2_train.faa",
        "dataset_splits/rep_split2_val.faa",
        "dataset_splits/rep_split2_test.faa",
    ],
    "split1_20pct": [
        "dataset_splits/rep_split1_20pct_train.faa",
        "dataset_splits/rep_split1_20pct_val.faa",
        "dataset_splits/rep_split1_20pct_test.faa",
    ],
    "split2_20pct": [
        "dataset_splits/rep_split2_20pct_train.faa",
        "dataset_splits/rep_split2_20pct_val.faa",
        "dataset_splits/rep_split2_20pct_test.faa",
    ],
}


def parse_faa(filepath):
    seqs = {}
    current_id  = None
    current_seq = []
    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_id is not None:
                    seqs[current_id] = ''.join(current_seq)
                parts       = line[1:].split('||')
                current_id  = parts[1].split()[0].strip()
                current_seq = []
            else:
                current_seq.append(line)
    if current_id is not None:
        seqs[current_id] = ''.join(current_seq)
    return seqs


def load_already_embedded(output_path):
    try:
        with h5py.File(output_path, 'r') as f:
            return set(f.keys())
    except Exception:
        return set()


def append_to_h5(output_path, batch_result):
    batch_dict = batch_result.to_dict()
    with h5py.File(output_path, 'a') as f:
        for seq_id, embedding in batch_dict.items():
            if seq_id not in f:
                f.create_dataset(seq_id, data=np.array(embedding))


def _try_embed(api, batch_dict, max_retries=3):
    for attempt in range(1, max_retries + 1):
        try:
            return api.embed(
                embedder_name=ESM2_650M,
                reduce=True,
                sequence_data=batch_dict,
                use_half_precision=False,
            ).run_with_progress()
        except Exception as e:
            if attempt == max_retries:
                raise
            wait = 2 ** attempt
            print(f"    Failed (attempt {attempt}/{max_retries}): {e} — retrying in {wait}s ...")
            time.sleep(wait)


def embed_batch_with_fallback(api, batch_dict):
    try:
        return _try_embed(api, batch_dict, max_retries=5)
    except Exception as e:
        print(f"  Full batch of {len(batch_dict)} failed: {e}")
        print(f"  Falling back to sub-batches of 100 ...")

    combined = None
    for sub_batch in batched(list(batch_dict.items()), 100):
        sub_dict = dict(sub_batch)
        try:
            result = _try_embed(api, sub_dict, max_retries=3)
        except Exception as e:
            print(f"    Sub-batch of 100 failed: {e} — falling back to groups of 10 ...")
            result = None

            for mini_batch in batched(sub_batch, 10):
                mini_dict = dict(mini_batch)
                try:
                    mini_result = _try_embed(api, mini_dict, max_retries=3)
                    result = mini_result if result is None else result.merge(mini_result)
                except Exception as e2:
                    print(f"    Skipping {len(mini_dict)} sequences: {list(mini_dict.keys())}")

        if result is not None:
            combined = result if combined is None else combined.merge(result)
        time.sleep(5)

    return combined


def embed_split(api, split_name, faa_files):
    output_path = os.path.join(OUT_DIR, f"{split_name}_esm2.h5")

    all_seqs = {}
    for faa in faa_files:
        all_seqs.update(parse_faa(faa))

    # Skip already-embedded sequences (resume support)
    already_done = load_already_embedded(output_path)
    if already_done:
        print(f"  Resuming — {len(already_done)} already saved, skipping those.")
        all_seqs = {k: v for k, v in all_seqs.items() if k not in already_done}

    if not all_seqs:
        print(f"  All sequences already embedded — skipping.")
        return

    seq_items    = list(all_seqs.items())
    total        = len(seq_items)
    n_batches    = (total + BATCH_SIZE - 1) // BATCH_SIZE
    print(f"  {total} sequences to embed ({n_batches} batches) ...")

    os.makedirs(OUT_DIR, exist_ok=True)
    for batch_idx, batch in enumerate(batched(seq_items, BATCH_SIZE), start=1):
        batch_dict = dict(batch)
        print(f"  Batch {batch_idx}/{n_batches} ({len(batch_dict)} seqs) ...", end=" ", flush=True)

        t0     = time.time()
        result = embed_batch_with_fallback(api, batch_dict)
        elapsed = time.time() - t0

        if result is not None:
            append_to_h5(output_path, result)
            n_saved = len(load_already_embedded(output_path))
            print(f"done ({elapsed:.1f}s)  [saved so far: {n_saved}]")
        else:
            print(f"batch returned no results — skipped.")

        time.sleep(10)  # rate limit buffer between batches

    print(f"  Saved to {output_path}")


if __name__ == "__main__":

    print("Connecting to biocentral API ...")
    api = BiocentralAPI().wait_until_healthy()

    for split_name, faa_files in SPLITS.items():
        print(f"\n{split_name} [ESM2 650M]")
        embed_split(api, split_name, faa_files)

    print("\nDone.")
