# ESM2 650M embeddings for the filtered non-extremophile sequences (baseline for the binary classifier), same
# procedure as compute_nonextremophile_embeddings.py. Duplicate IDs were added later by
# compute_baseline_embeddings.py.

import time
import h5py
import numpy as np
from pathlib import Path
from biocentral_api import BiocentralAPI, CommonEmbedder, batched

BASE_DIR  = Path(__file__).resolve().parents[2]
INPUT_FAA = BASE_DIR / "nonextremophile_data" / "nonextremophile_filtered.faa"
OUT_H5    = BASE_DIR / "embeddings" / "nonextremophile_esm2.h5"

ESM2_650M = CommonEmbedder("facebook/esm2_t33_650M_UR50D")


def parse_faa(filepath):
    raw = {}
    current_id  = None
    current_seq = []

    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if current_id is not None:
                    raw[current_id] = "".join(current_seq)
                current_id  = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)

    if current_id is not None:
        raw[current_id] = "".join(current_seq)

    seen_seqs = {}
    for pid, seq in raw.items():
        if seq not in seen_seqs:
            seen_seqs[seq] = pid

    sequences = {pid: seq for seq, pid in seen_seqs.items()}
    n_dupes = len(raw) - len(sequences)
    if n_dupes:
        print(f"  Removed {n_dupes:,} duplicate sequences.")

    return sequences


def load_already_embedded(output_path):
    try:
        with h5py.File(output_path, "r") as f:
            return set(f.keys())
    except Exception:
        return set()


def append_to_h5(output_path, batch_result):
    batch_dict = batch_result.to_dict()
    with h5py.File(output_path, "a") as f:
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
            print(f"    Attempt {attempt}/{max_retries} failed: {e} — retrying in {wait}s ...")
            time.sleep(wait)


def embed_batch_with_retry(api, batch_dict):
    try:
        return _try_embed(api, batch_dict, max_retries=5)
    except Exception as e:
        print(f"  Full batch of {len(batch_dict)} failed: {e} — falling back to sub-batches ...")

    items    = list(batch_dict.items())
    combined = None
    for sub_batch in batched(items, 100):
        sub_dict = dict(sub_batch)
        try:
            result = _try_embed(api, sub_dict, max_retries=3)
        except Exception as e:
            print(f"    Sub-batch failed: {e} — trying groups of 10 ...")
            result = None
            for mini_batch in batched(sub_batch, 10):
                mini_dict = dict(mini_batch)
                try:
                    mini_result = _try_embed(api, mini_dict, max_retries=3)
                    result = mini_result if result is None else result.merge(mini_result)
                except Exception as e2:
                    print(f"    Skipping {len(mini_dict)} sequences permanently.")
        if result is not None:
            combined = result if combined is None else combined.merge(result)
        time.sleep(5)
    return combined


if __name__ == "__main__":
    OUT_H5.parent.mkdir(exist_ok=True)

    print(f"Parsing {INPUT_FAA} ...")
    all_seqs = parse_faa(INPUT_FAA)
    print(f"Total unique sequences: {len(all_seqs):,}")

    already_done = load_already_embedded(OUT_H5)
    if already_done:
        print(f"Resuming — {len(already_done):,} already embedded, skipping.")
        all_seqs = {k: v for k, v in all_seqs.items() if k not in already_done}
        print(f"Remaining: {len(all_seqs):,}")

    if not all_seqs:
        print("All sequences already embedded.")
    else:
        print("\nConnecting to biocentral API ...")
        api = BiocentralAPI().wait_until_healthy()

        seq_items     = list(all_seqs.items())
        total_batches = (len(seq_items) + 499) // 500

        for batch_idx, batch in enumerate(batched(seq_items, 500), start=1):
            print(f"\nEmbedding batch {batch_idx}/{total_batches} ...")
            result = embed_batch_with_retry(api, dict(batch))
            if result is not None:
                append_to_h5(OUT_H5, result)
                print(f"  Saved batch {batch_idx} to {OUT_H5}")
            time.sleep(10)

    print(f"\nDone. Embeddings saved to {OUT_H5}")
