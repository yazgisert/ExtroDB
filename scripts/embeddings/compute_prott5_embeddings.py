# ProtT5 embeddings for the 5-class (split1) and 12-class (split2) splits at 30% and 20% identity via the
# biocentral API. Batches of 500, falling back to 100 and then 10 if a batch fails; resumes from the
# existing .h5 files.
# Run: caffeinate -i python scripts/embeddings/compute_prott5_embeddings.py

import time
import h5py
import numpy as np
from biocentral_api import BiocentralAPI, CommonEmbedder, batched


def parse_faa(filepath):
    # Key is the protein accession; the proteome ID is shared by many proteins.
    sequences = {}
    current_id = None
    current_seq = []

    with open(filepath) as f:
        for line in f:
            line = line.strip()

            if line.startswith('>'):
                if current_id is not None:
                    sequences[current_id] = ''.join(current_seq)

                # Header format: >proteome_id||protein_id description||label
                parts = line[1:].split('||')
                current_id = parts[1].split()[0].strip()
                current_seq = []

            else:
                current_seq.append(line)

    if current_id is not None:
        sequences[current_id] = ''.join(current_seq)

    return sequences


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
                embedder_name=CommonEmbedder.ProtT5,
                reduce=True,
                sequence_data=batch_dict,
                use_half_precision=False,
            ).run_with_progress()
        except Exception as e:
            if attempt == max_retries:
                raise
            wait = 2 ** attempt  # 2s, 4s, 8s ...
            print(f"    Failed (attempt {attempt}/{max_retries}): {e} — retrying in {wait}s ...")
            time.sleep(wait)


def embed_batch_with_retry(api, batch_dict):
    try:
        return _try_embed(api, batch_dict, max_retries=5)
    except Exception as e:
        print(f"  Full batch of {len(batch_dict)} failed permanently: {e}")
        print(f"  Falling back to sub-batches of 100 ...")

    items = list(batch_dict.items())
    combined = None
    for sub_batch in batched(items, 100):
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
                    print(f"    Skipping {len(mini_dict)} sequences permanently: {list(mini_dict.keys())}")

        if result is not None:
            combined = result if combined is None else combined.merge(result)
        time.sleep(5)

    return combined


def embed_faa_files(faa_files, output_path):

    all_seqs = {}
    for faa_file in faa_files:
        print(f"Parsing {faa_file} ...")
        parsed = parse_faa(faa_file)
        all_seqs.update(parsed)
    print(f"Total unique sequences: {len(all_seqs)}")

    already_done = load_already_embedded(output_path)
    if already_done:
        print(f"Resuming — {len(already_done)} sequences already saved, skipping those.")
        all_seqs = {k: v for k, v in all_seqs.items() if k not in already_done}
        print(f"Remaining sequences to embed: {len(all_seqs)}")

    if not all_seqs:
        print("All sequences already embedded. Nothing to do.")
        return

    print("Connecting to biocentral API ...")
    api = BiocentralAPI().wait_until_healthy()

    seq_items = list(all_seqs.items())
    total_batches = (len(seq_items) + 499) // 500

    for batch_idx, batch in enumerate(batched(seq_items, 500), start=1):
        print(f"Embedding batch {batch_idx}/{total_batches} ...")

        batch_dict = {seq_id: seq for seq_id, seq in batch}
        batch_result = embed_batch_with_retry(api, batch_dict)

        if batch_result is not None:
            append_to_h5(output_path, batch_result)
            print(f"  Saved batch {batch_idx} to {output_path}")

        # Pause between batches to stay under rate limits
        time.sleep(10)

    print(f"\nDone. All embeddings saved to {output_path}")


if __name__ == "__main__":

    import os
    os.makedirs("embeddings", exist_ok=True)

    # One .h5 per split (train + val + test together).

    runs = [
        (["dataset_splits/rep_split1_train.faa",
          "dataset_splits/rep_split1_val.faa",
          "dataset_splits/rep_split1_test.faa"],       "embeddings/split1_30pct.h5"),

        (["dataset_splits/rep_split2_train.faa",
          "dataset_splits/rep_split2_val.faa",
          "dataset_splits/rep_split2_test.faa"],       "embeddings/split2_30pct.h5"),

        (["dataset_splits/rep_split1_20pct_train.faa",
          "dataset_splits/rep_split1_20pct_val.faa",
          "dataset_splits/rep_split1_20pct_test.faa"], "embeddings/split1_20pct.h5"),

        (["dataset_splits/rep_split2_20pct_train.faa",
          "dataset_splits/rep_split2_20pct_val.faa",
          "dataset_splits/rep_split2_20pct_test.faa"], "embeddings/split2_20pct.h5"),
    ]

    for faa_files, output_path in runs:
        print(f"\n{output_path}")
        embed_faa_files(faa_files, output_path)
