# Lists the leakage-filtered algae proteins whose sequence is identical to another algae protein. One protein per
# sequence is kept, preferably one with an embedding, then by ID; the others are written out and skipped by the algae
# evaluations.

import h5py

from eval_eukaryote_embedders import ALGAE_EMBEDDINGS, ALGAE_FAA, BLACKLIST, read_faa

OUT_TXT = "results/eukaryote_leakage/duplicate_protein_ids.txt"


if __name__ == "__main__":
    with open(BLACKLIST) as f:
        blacklist = {line.strip() for line in f if line.strip()}
    with h5py.File(ALGAE_EMBEDDINGS["prott5"], "r") as f:
        embedded = set(f.keys())
    seqs = read_faa(ALGAE_FAA)

    seen, duplicates = set(), []
    for protein_id in sorted(seqs, key=lambda p: (p not in embedded, p)):
        if protein_id in blacklist:
            continue
        if seqs[protein_id] in seen:
            duplicates.append(protein_id)
        seen.add(seqs[protein_id])

    with open(OUT_TXT, "w") as f:
        f.write("\n".join(duplicates) + "\n")
    print(f"{len(duplicates):,} duplicate proteins of {len(seen):,} unique sequences")
    print(f"Saved to {OUT_TXT}")
