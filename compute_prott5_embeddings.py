# compute_prott5_embeddings.py
#
# Computes ProtT5 embeddings for the representative sequence splits
# using the biocentral API and the ProtT5 protein language model.
#
import time
import h5py
import numpy as np
from biocentral_api import BiocentralAPI, CommonEmbedder, batched


def parse_faa(filepath):
    """
    Reads one of our custom .faa files and returns a plain dictionary:
        { sequence_id: amino_acid_sequence }

    Our FASTA headers look like:
        >GCF_003151135.1||WP_008296005.1 MULTISPECIES: hypothetical protein||halophile

    We use the protein accession (e.g. WP_008296005.1) as the key because
    multiple proteins share the same proteome ID — using that would overwrite entries.
    """
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
                current_id = parts[1].split()[0].strip()  # first word of second segment
                current_seq = []

            else:
                current_seq.append(line)

    if current_id is not None:
        sequences[current_id] = ''.join(current_seq)

    return sequences


def load_already_embedded(output_path):
    """
    If a checkpoint .h5 file already exists, return the set of sequence IDs
    already saved so we can skip them on resume.
    """
    try:
        with h5py.File(output_path, 'r') as f:
            return set(f.keys())
    except Exception:
        return set()


def append_to_h5(output_path, batch_result):
    """
    Append a batch's embeddings directly to the .h5 file using h5py in append
    mode. This means we never overwrite previously saved sequences — each restart
    safely adds on top of what's already there.
    """
    batch_dict = batch_result.to_dict()  # {seq_id: np.array}
    with h5py.File(output_path, 'a') as f:  # 'a' = append/create mode
        for seq_id, embedding in batch_dict.items():
            if seq_id not in f:  # skip if somehow already present
                f.create_dataset(seq_id, data=np.array(embedding))


def _try_embed(api, batch_dict, max_retries=3):
    """
    Try to embed a dict of sequences, retrying on failure with exponential backoff.
    Returns the result, or raises if all retries are exhausted.
    """
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
    """
    Embed a batch with automatic fallback:
      1. Try the full batch (retries on transient errors).
      2. If it keeps failing, split into sub-batches of 100 and try each.
      3. If a sub-batch of 100 fails, split into groups of 10 and try each.
      4. If a group of 10 still fails, log and skip those sequences.
    """
    # 1. Try the full batch
    try:
        return _try_embed(api, batch_dict, max_retries=5)
    except Exception as e:
        print(f"  Full batch of {len(batch_dict)} failed permanently: {e}")
        print(f"  Falling back to sub-batches of 100 ...")

    # 2. Fall back to sub-batches of 100
    items = list(batch_dict.items())
    combined = None
    for sub_batch in batched(items, 100):
        sub_dict = dict(sub_batch)
        try:
            result = _try_embed(api, sub_dict, max_retries=3)
        except Exception as e:
            print(f"    Sub-batch of 100 failed: {e} — falling back to groups of 10 ...")
            result = None

            # 3. Fall back to groups of 10
            for mini_batch in batched(sub_batch, 10):
                mini_dict = dict(mini_batch)
                try:
                    mini_result = _try_embed(api, mini_dict, max_retries=3)
                    result = mini_result if result is None else result.merge(mini_result)
                except Exception as e2:
                    # 4. Give up, log skipped IDs
                    print(f"    Skipping {len(mini_dict)} sequences permanently: {list(mini_dict.keys())}")

        if result is not None:
            combined = result if combined is None else combined.merge(result)
        time.sleep(5)  # pause between sub-batches

    return combined


def embed_faa_files(faa_files, output_path):
    """
    Reads one or more .faa files, embeds all sequences using ProtT5 via the
    biocentral API, and saves results to a single .h5 file.

    Each batch is appended to the .h5 file immediately after it succeeds, so
    progress is never lost on a crash or restart. Resume is automatic.
    """

    # Step 1: Parse all files into one combined {protein_id: sequence} dict
    all_seqs = {}
    for faa_file in faa_files:
        print(f"Parsing {faa_file} ...")
        parsed = parse_faa(faa_file)
        all_seqs.update(parsed)
    print(f"Total unique sequences: {len(all_seqs)}")

    # Step 2: Skip sequences already saved in a previous run
    already_done = load_already_embedded(output_path)
    if already_done:
        print(f"Resuming — {len(already_done)} sequences already saved, skipping those.")
        all_seqs = {k: v for k, v in all_seqs.items() if k not in already_done}
        print(f"Remaining sequences to embed: {len(all_seqs)}")

    if not all_seqs:
        print("All sequences already embedded. Nothing to do.")
        return

    # Step 3: Connect to the biocentral API
    print("Connecting to biocentral API ...")
    api = BiocentralAPI().wait_until_healthy()

    # Step 4: Embed in batches of 500 and append each batch to the .h5 file immediately
    seq_items = list(all_seqs.items())
    total_batches = (len(seq_items) + 499) // 500

    for batch_idx, batch in enumerate(batched(seq_items, 500), start=1):
        print(f"Embedding batch {batch_idx}/{total_batches} ...")

        batch_dict = {seq_id: seq for seq_id, seq in batch}
        batch_result = embed_batch_with_retry(api, batch_dict)

        if batch_result is not None:
            # Append this batch directly to the .h5 file — never overwrites existing data
            append_to_h5(output_path, batch_result)
            print(f"  Saved batch {batch_idx} to {output_path}")

        # Pause between batches to stay under rate limits
        time.sleep(10)

    print(f"\nDone. All embeddings saved to {output_path}")


if __name__ == "__main__":

    import os
    os.makedirs("embeddings", exist_ok=True)

    # Each split (30% and 20% clustering threshold) has two label sets:
    #   split1 = 5 coarse classes, split2 = 12 fine-grained classes.
    # We embed train + val + test together into one .h5 file per split.

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
        print(f"\n{'='*60}")
        print(f"Starting: {output_path}")
        print(f"{'='*60}")
        embed_faa_files(faa_files, output_path)
