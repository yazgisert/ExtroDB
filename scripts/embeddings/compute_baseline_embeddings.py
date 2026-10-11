# ESM2 650M / ESM2 8M embeddings for the baseline comparison via the biocentral API: 5- and 12-class splits,
# non-extremophiles and leakage-filtered algae (same proteins as the ProtT5 files). Identical sequences are
# embedded once and stored under every ID.
# Run: python scripts/embeddings/compute_baseline_embeddings.py --embedder esm2_8m --dataset splits

import argparse
import os
import time
from collections import defaultdict

import h5py
import numpy as np
from biocentral_api import BiocentralAPI, CommonEmbedder, batched

BATCH_SIZE = 1000

EMBEDDERS = {
    "esm2_650m": (CommonEmbedder.ESM2_650M, "_esm2"),
    "esm2_8m":   (CommonEmbedder.ESM2_8M,   "_esm2_8m"),
}

SPLITS = {
    "split1_30pct": ["dataset_splits/rep_split1_train.faa",
                     "dataset_splits/rep_split1_val.faa",
                     "dataset_splits/rep_split1_test.faa"],
    "split2_30pct": ["dataset_splits/rep_split2_train.faa",
                     "dataset_splits/rep_split2_val.faa",
                     "dataset_splits/rep_split2_test.faa"],
    "split1_20pct": ["dataset_splits/rep_split1_20pct_train.faa",
                     "dataset_splits/rep_split1_20pct_val.faa",
                     "dataset_splits/rep_split1_20pct_test.faa"],
    "split2_20pct": ["dataset_splits/rep_split2_20pct_train.faa",
                     "dataset_splits/rep_split2_20pct_val.faa",
                     "dataset_splits/rep_split2_20pct_test.faa"],
}

NONEXTREM_FAA = "nonextremophile_data/nonextremophile_filtered.faa"
ALGAE_FAA     = "results/homology_inference/filtered_eukaryote_query.faa"
BLACKLIST     = "results/eukaryote_leakage/flagged_protein_ids.txt"


def parse_faa(filepath, id_from_header):
    seqs = {}
    current_id, current_seq = None, []
    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if current_id is not None:
                    seqs[current_id] = "".join(current_seq)
                current_id, current_seq = id_from_header(line[1:]), []
            else:
                current_seq.append(line)
    if current_id is not None:
        seqs[current_id] = "".join(current_seq)
    return seqs


def split_id(header):
    return header.split("||")[1].split()[0].strip()


def first_word(header):
    return header.split()[0]


def h5_keys(path):
    try:
        with h5py.File(path, "r") as f:
            return set(f.keys())
    except (FileNotFoundError, OSError):
        return set()


def embed_with_retry(api, embedder, batch, max_retries=3):
    for attempt in range(1, max_retries + 1):
        try:
            return api.embed(
                embedder_name=embedder,
                reduce=True,
                sequence_data=batch,
                use_half_precision=False,
            ).run_with_progress()
        except Exception as e:
            if attempt == max_retries:
                raise
            print(f"    attempt {attempt}/{max_retries} failed: {e} — retrying in {2 ** attempt}s")
            time.sleep(2 ** attempt)


def embed_batch(api, embedder, batch):
    try:
        return embed_with_retry(api, embedder, batch, max_retries=5)
    except Exception as e:
        print(f"  batch of {len(batch)} failed ({e}), trying smaller batches ...")

    combined = None
    for size in (100, 10):
        failed = {}
        for sub in batched(list(batch.items()), size):
            sub = dict(sub)
            try:
                result = embed_with_retry(api, embedder, sub)
                combined = result if combined is None else combined.merge(result)
            except Exception:
                failed.update(sub)
            time.sleep(2)
        if not failed:
            break
        batch = failed
    else:
        print(f"    skipping {len(failed)} sequences: {list(failed)[:5]} ...")
    return combined


def save(path, result, ids_by_seq, seqs):
    with h5py.File(path, "a") as f:
        for rep_id, embedding in result.to_dict().items():
            for seq_id in ids_by_seq[seqs[rep_id]]:
                if seq_id not in f:
                    f.create_dataset(seq_id, data=np.array(embedding))


def embed_all(api, embedder, seqs, out_path, earlier=()):
    # The API rejects duplicate sequences, so each unique sequence is embedded once
    # and stored under every protein ID that shares it. Sequences already embedded
    # in an earlier file of the same run (e.g. another split) are copied over.
    ids_by_seq = defaultdict(list)
    for seq_id, seq in seqs.items():
        ids_by_seq[seq].append(seq_id)

    elsewhere = {}
    for other_seqs, other_path in earlier:
        other_done = h5_keys(other_path)
        for seq_id, seq in other_seqs.items():
            if seq_id in other_done:
                elsewhere.setdefault(seq, (other_path, seq_id))

    done = h5_keys(out_path)
    copied = 0
    todo = {}
    others = {path: h5py.File(path, "r") for _, path in earlier}
    with h5py.File(out_path, "a") as f:
        for seq, ids in ids_by_seq.items():
            have = [i for i in ids if i in done]
            if len(have) == len(ids):
                continue
            if have:
                embedding = f[have[0]][:]
            elif seq in elsewhere:
                other_path, other_id = elsewhere[seq]
                embedding = others[other_path][other_id][:]
            else:
                todo[ids[0]] = seq
                continue
            for i in ids:
                if i not in f:
                    f.create_dataset(i, data=embedding)
                    copied += 1
    for other in others.values():
        other.close()

    print(f"  {len(seqs):,} proteins ({len(ids_by_seq):,} unique sequences), "
          f"{len(done & seqs.keys()):,} already embedded, {copied:,} copied from identical sequences, "
          f"{len(todo):,} unique sequences to embed")
    if not todo:
        return

    n_batches = (len(todo) + BATCH_SIZE - 1) // BATCH_SIZE
    for i, batch in enumerate(batched(list(todo.items()), BATCH_SIZE), start=1):
        t0 = time.time()
        result = embed_batch(api, embedder, dict(batch))
        if result is not None:
            save(out_path, result, ids_by_seq, seqs)
        print(f"  batch {i}/{n_batches} done ({time.time() - t0:.0f}s)")
        time.sleep(2)
    print(f"  Saved to {out_path}")


def jobs_for(dataset, suffix):
    if dataset == "splits":
        for name, files in SPLITS.items():
            seqs = {}
            for faa in files:
                seqs.update(parse_faa(faa, split_id))
            yield name, seqs, f"embeddings/{name}{suffix}.h5"

    elif dataset == "nonextremophile":
        # Only the proteins in the ProtT5 file, so every embedder uses the same set
        prott5_ids = h5_keys("embeddings/nonextremophile.h5")
        seqs = parse_faa(NONEXTREM_FAA, first_word)
        yield dataset, {k: v for k, v in seqs.items() if k in prott5_ids}, f"embeddings/nonextremophile{suffix}.h5"

    elif dataset == "algae":
        prott5_ids = h5_keys("embeddings/eukaryote_algae.h5")
        with open(BLACKLIST) as f:
            blacklist = {line.strip() for line in f if line.strip()}
        seqs = parse_faa(ALGAE_FAA, first_word)
        keep = {k: v for k, v in seqs.items() if k in prott5_ids and k not in blacklist}
        yield dataset, keep, f"embeddings/eukaryote_algae_filtered{suffix}.h5"


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--embedder", choices=EMBEDDERS, required=True)
    parser.add_argument("--dataset", choices=["splits", "nonextremophile", "algae"], required=True)
    args = parser.parse_args()

    embedder, suffix = EMBEDDERS[args.embedder]
    os.makedirs("embeddings", exist_ok=True)

    print("Connecting to biocentral API ...")
    api = BiocentralAPI().wait_until_healthy()

    earlier = []
    for name, seqs, out_path in jobs_for(args.dataset, suffix):
        print(f"\n{name} [{args.embedder}]")
        embed_all(api, embedder, seqs, out_path, earlier)
        earlier.append((seqs, out_path))

    print("\nDone.")
