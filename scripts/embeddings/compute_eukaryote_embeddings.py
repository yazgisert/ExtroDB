# ProtT5 embeddings for the 13 algae proteomes via the biocentral API, plus the label table (protein_id,
# proteome_id, broad_class). Identical sequences are embedded once; batches of 500 fall back to 100 and then
# 10 if they fail; resumes from the .h5 file.

import csv
import time
import h5py
import numpy as np
from pathlib import Path
from biocentral_api import BiocentralAPI, CommonEmbedder, batched

BASE_DIR    = Path(__file__).resolve().parents[2]
INDEX_PATH  = BASE_DIR / "searchable_resource" / "extremophile_master_index.tsv"
OUT_H5      = BASE_DIR / "embeddings" / "eukaryote_algae.h5"
OUT_LABELS  = BASE_DIR / "embeddings" / "eukaryote_algae_labels.tsv"

CLASS_MAP = {
    "Thermophile":            "thermophile",
    "Hyperthermophile":       "thermophile",
    "Thermophilic Alga":      "thermophile",
    "High-Temperature Alga":  "thermophile",
    "Psychrophile":           "psychrophile",
    "Psychrophilic Alga":     "psychrophile",
    "Halophile (Slight)":     "halophile",
    "Halophile (Moderate)":   "halophile",
    "Halophile (Extreme)":    "halophile",
    "Halophilic Alga":        "halophile",
    "Saline-Adapted Alga":    "halophile",
    "Acidophile":             "acidophile",
    "Moderate Acidophile":    "acidophile",
    "Hyperacidophile":        "acidophile",
    "Acidotolerant":          "acidophile",
    "Acidophilic Alga":       "acidophile",
    "Thermoacidophilic Alga": "acidophile",   # primary = acidophile for multi-label
    "Alkaliphile":            "alkaliphile",
    "Extreme Alkaliphile":    "alkaliphile",
    "Alkalitolerant":         "alkaliphile",
}

CLASS_ORDER = ["thermophile", "halophile", "acidophile", "psychrophile", "alkaliphile"]


def get_broad_class(phenotype_labels_str, phenotype_focus):
    # phenotype_focus decides for multi-label organisms, otherwise the first mappable label.
    if phenotype_focus and phenotype_focus in CLASS_ORDER:
        return phenotype_focus

    labels = [l.strip() for l in phenotype_labels_str.split("|")]
    for label in labels:
        if label in CLASS_MAP:
            return CLASS_MAP[label]
    return None


def load_algae_index():
    with open(INDEX_PATH) as f:
        rows = list(csv.DictReader(f, delimiter="\t"))

    algae = []
    for row in rows:
        if row.get("domain_scope") != "eukaryote":
            continue
        broad = get_broad_class(row["phenotype_labels"], row.get("phenotype_focus", ""))
        if broad is None:
            print(f"  WARNING: no class mapping for {row['taxon_name']} — skipping")
            continue
        fasta_path_str = (row.get("local_fasta_path") or "").strip()
        if not fasta_path_str:
            print(f"  WARNING: no FASTA path recorded for {row['taxon_name']} — skipping")
            continue
        fasta_path = BASE_DIR / fasta_path_str
        if not fasta_path.exists():
            print(f"  WARNING: FASTA not found: {fasta_path} — skipping")
            continue
        algae.append({
            "proteome_id": row["proteome_id"],
            "taxon_name":  row["taxon_name"],
            "broad_class": broad,
            "fasta_path":  fasta_path,
        })

    print(f"Loaded {len(algae)} eukaryote organisms from master index.")
    return algae


def parse_ncbi_faa(fasta_path, proteome_id, broad_class):
    sequences = {}
    meta      = {}
    current_id  = None
    current_seq = []

    with open(fasta_path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if current_id is not None:
                    sequences[current_id] = "".join(current_seq)
                    meta[current_id] = (proteome_id, broad_class)
                current_id  = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)

    if current_id is not None:
        sequences[current_id] = "".join(current_seq)
        meta[current_id] = (proteome_id, broad_class)

    return sequences, meta


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
                embedder_name=CommonEmbedder.ProtT5,
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

    items = list(batch_dict.items())
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


def main():
    OUT_H5.parent.mkdir(exist_ok=True)

    algae = load_algae_index()

    all_seqs = {}
    all_meta = {}
    for org in algae:
        print(f"  Parsing {org['taxon_name']} ({org['broad_class']}) ...")
        seqs, meta = parse_ncbi_faa(org["fasta_path"], org["proteome_id"], org["broad_class"])
        all_seqs.update(seqs)
        all_meta.update(meta)

    print(f"\nTotal sequences: {len(all_seqs)}")

    # Identical sequences give identical embeddings, so each is embedded once.
    seen_seqs: dict[str, str] = {}
    for pid, seq in all_seqs.items():
        if seq not in seen_seqs:
            seen_seqs[seq] = pid
    deduped = {pid: seq for seq, pid in seen_seqs.items()}
    n_dupes = len(all_seqs) - len(deduped)
    if n_dupes:
        print(f"  Removed {n_dupes:,} duplicate sequences (same amino acid content across organisms).")
    all_seqs = deduped

    already_done = load_already_embedded(OUT_H5)
    if already_done:
        print(f"Resuming — {len(already_done)} already embedded, skipping.")
        all_seqs = {k: v for k, v in all_seqs.items() if k not in already_done}
        print(f"Remaining: {len(all_seqs)}")

    if not all_seqs:
        print("All sequences already embedded.")
    else:
        print("\nConnecting to biocentral API ...")
        api = BiocentralAPI().wait_until_healthy()

        seq_items     = list(all_seqs.items())
        total_batches = (len(seq_items) + 499) // 500

        for batch_idx, batch in enumerate(batched(seq_items, 500), start=1):
            print(f"\nEmbedding batch {batch_idx}/{total_batches} ...")
            batch_dict = dict(batch)
            result = embed_batch_with_retry(api, batch_dict)
            if result is not None:
                append_to_h5(OUT_H5, result)
                print(f"  Saved batch {batch_idx} to {OUT_H5}")
            time.sleep(10)

    print(f"\nSaving label mapping to {OUT_LABELS}")
    with open(OUT_LABELS, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["protein_id", "proteome_id", "broad_class"])
        for protein_id, (proteome_id, broad_class) in all_meta.items():
            writer.writerow([protein_id, proteome_id, broad_class])

    print(f"\nDone. Embeddings: {OUT_H5}  |  Labels: {OUT_LABELS}")


if __name__ == "__main__":
    main()
