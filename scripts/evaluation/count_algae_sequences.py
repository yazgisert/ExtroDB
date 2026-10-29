# Counts the algae test set per species by unique sequences: proteins, unique sequences, unique sequences removed by
# the leakage filter (check_eukaryote_leakage.py) and unique sequences used. Identical sequences have identical hits,
# so a sequence is removed if any of its proteins was flagged.

import csv
import json

MASTER_INDEX = "searchable_resource/extremophile_master_index.tsv"
FLAGGED_IDS  = "results/eukaryote_leakage/flagged_protein_ids.txt"
OUT_JSON     = "results/eukaryote_leakage/algae_unique_counts.json"


def read_faa(path):
    seqs, current = {}, None
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                current = line[1:].split()[0]
                seqs[current] = []
            else:
                seqs[current].append(line.strip())
    return {protein_id: "".join(chunks) for protein_id, chunks in seqs.items()}


if __name__ == "__main__":
    flagged = {line.strip() for line in open(FLAGGED_IDS) if line.strip()}
    algae = [r for r in csv.DictReader(open(MASTER_INDEX), delimiter="\t")
             if r["domain_scope"] == "eukaryote" and r["proteome_status"] == "downloaded"]

    result, all_seqs, removed_seqs = {"species": {}}, set(), set()
    print(f"{'species':<32} {'proteins':>9} {'unique':>8} {'removed':>8} {'used':>8}")
    for r in sorted(algae, key=lambda r: r["taxon_name"]):
        seqs = read_faa(r["local_fasta_path"])
        removed = {seq for protein_id, seq in seqs.items() if protein_id in flagged}
        unique = set(seqs.values())
        all_seqs |= unique
        removed_seqs |= removed
        result["species"][r["taxon_name"]] = {"proteins": len(seqs), "unique": len(unique), "removed": len(removed),
                                              "used": len(unique - removed)}
        print(f"{r['taxon_name']:<32} {len(seqs):>9,} {len(unique):>8,} {len(removed):>8,} {len(unique - removed):>8,}")

    result["total"] = {"proteins": sum(s["proteins"] for s in result["species"].values()), "unique": len(all_seqs),
                       "removed": len(removed_seqs), "used": len(all_seqs - removed_seqs)}
    t = result["total"]
    print(f"{'total (distinct sequences)':<32} {t['proteins']:>9,} {t['unique']:>8,} {t['removed']:>8,} {t['used']:>8,}")

    with open(OUT_JSON, "w") as f:
        json.dump(result, f, indent=2)
    print(f"\nSaved to {OUT_JSON}")
