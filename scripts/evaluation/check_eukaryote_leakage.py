# Flags algae proteins with >= 30% identity (MMseqs2) to any prokaryote train/val protein; the flagged IDs
# are left out of all eukaryote evaluations.

import csv
import subprocess
import tempfile
from collections import defaultdict
from pathlib import Path


BASE_DIR   = Path(__file__).resolve().parents[2]
SPLITS_DIR = BASE_DIR / "dataset_splits"
INDEX_PATH = BASE_DIR / "searchable_resource" / "extremophile_master_index.tsv"

PROKARYOTE_REFS = [
    SPLITS_DIR / "rep_split1_train.faa",
    SPLITS_DIR / "rep_split1_val.faa",
]

OUT_DIR      = BASE_DIR / "results" / "eukaryote_leakage"
OUT_HITS     = OUT_DIR / "leakage_hits.tsv"
OUT_BLACKLIST = OUT_DIR / "flagged_protein_ids.txt"
OUT_SUMMARY  = OUT_DIR / "leakage_summary.txt"

IDENTITY_THRESHOLD = 0.30
SENSITIVITY        = 7.5


def load_eukaryote_organisms():
    organisms = []
    with open(INDEX_PATH) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            if row.get("domain_scope") != "eukaryote":
                continue
            fasta_path_str = (row.get("local_fasta_path") or "").strip()
            if not fasta_path_str:
                print(f"  WARNING: no FASTA path for {row['taxon_name']} — skipping")
                continue
            fasta_path = BASE_DIR / fasta_path_str
            if not fasta_path.exists():
                print(f"  WARNING: FASTA not found: {fasta_path} — skipping")
                continue
            organisms.append({
                "taxon_name":  row["taxon_name"],
                "proteome_id": row["proteome_id"],
                "fasta_path":  fasta_path,
                "broad_class": row.get("phenotype_focus") or row.get("phenotype_labels", "unknown"),
            })
    print(f"Found {len(organisms)} eukaryote organisms in master index.")
    return organisms


def build_combined_eukaryote_fasta(organisms, out_path):
    protein_class = {}
    total = 0
    with open(out_path, "w") as out_fh:
        for org in organisms:
            with open(org["fasta_path"]) as fh:
                current_id = None
                for line in fh:
                    if line.startswith(">"):
                        current_id = line[1:].split()[0]
                        protein_class[current_id] = org["broad_class"]
                        total += 1
                    out_fh.write(line)
    print(f"Combined {total:,} eukaryote proteins into query FASTA.")
    return protein_class


def build_combined_prokaryote_fasta(ref_paths, out_path):
    total = 0
    with open(out_path, "w") as out_fh:
        for path in ref_paths:
            with open(path) as fh:
                for line in fh:
                    out_fh.write(line)
                    if line.startswith(">"):
                        total += 1
    print(f"Combined {total:,} prokaryote sequences (train + val) into reference FASTA.")
    return total


def run_mmseqs(query, reference, results_file, tmp_dir):
    cmd = [
        "mmseqs", "easy-search",
        str(query), str(reference), str(results_file), str(tmp_dir),
        "--min-seq-id", str(IDENTITY_THRESHOLD),
        "-s",          str(SENSITIVITY),
        "--format-output", "query,target,pident,evalue,bits",
        "--threads",   "4",
    ]
    print(f"\nRunning MMseqs2 (threshold: {IDENTITY_THRESHOLD*100:.0f}% identity)...")
    print("  " + " ".join(cmd))
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(result.stderr[-3000:])
        raise RuntimeError("MMseqs2 failed — see stderr above.")
    for line in result.stderr.splitlines():
        if line.strip():
            print(f"  [mmseqs] {line}")


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    for ref in PROKARYOTE_REFS:
        if not ref.exists():
            raise SystemExit(f"Missing reference FAA: {ref}")

    organisms = load_eukaryote_organisms()

    with tempfile.TemporaryDirectory() as tmp_str:
        tmp_dir   = Path(tmp_str)
        query_faa = tmp_dir / "eukaryote_query.faa"
        ref_faa   = tmp_dir / "prokaryote_ref.faa"
        hits_raw  = tmp_dir / "mmseqs_hits.tsv"

        print("\nBuilding query FASTA (eukaryote proteins)...")
        protein_class = build_combined_eukaryote_fasta(organisms, query_faa)
        total_proteins = len(protein_class)

        print("\nBuilding reference FASTA (prokaryote train + val)...")
        build_combined_prokaryote_fasta(PROKARYOTE_REFS, ref_faa)

        mmseqs_tmp = tmp_dir / "mmseqs_tmp"
        mmseqs_tmp.mkdir()
        run_mmseqs(query_faa, ref_faa, hits_raw, mmseqs_tmp)

        flagged_ids = set()
        hits = []
        with open(hits_raw) as fh:
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 5:
                    continue
                query_id, target_id, pident, evalue, bits = parts[:5]
                flagged_ids.add(query_id)
                hits.append((query_id, target_id, float(pident), evalue, bits))

        print(f"\nMMseqs2 found {len(flagged_ids):,} eukaryote proteins with >= {IDENTITY_THRESHOLD*100:.0f}% identity to a prokaryote training protein.")
        print(f"That is {len(flagged_ids)/total_proteins*100:.2f}% of the {total_proteins:,} total eukaryote proteins.")

        with open(OUT_HITS, "w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t")
            w.writerow(["eukaryote_protein", "prokaryote_training_protein", "pident_%", "evalue", "bits"])
            for row in sorted(hits, key=lambda x: -x[2]):
                w.writerow([row[0], row[1], f"{row[2]:.1f}", row[3], row[4]])
        print(f"Hit table saved: {OUT_HITS}")

        with open(OUT_BLACKLIST, "w") as fh:
            for pid in sorted(flagged_ids):
                fh.write(pid + "\n")
        print(f"Blacklist saved: {OUT_BLACKLIST}")

        flagged_by_class = defaultdict(int)
        total_by_class   = defaultdict(int)
        for pid, cls in protein_class.items():
            total_by_class[cls] += 1
            if pid in flagged_ids:
                flagged_by_class[cls] += 1

        print("\nFlagged proteins per class:")
        for cls in sorted(total_by_class):
            n_total   = total_by_class[cls]
            n_flagged = flagged_by_class[cls]
            print(f"  {cls:<20} {n_flagged:>6,} / {n_total:>6,}  ({n_flagged/n_total*100:.2f}%)")

        with open(OUT_SUMMARY, "w") as fh:
            fh.write("Eukaryote Leakage Check Summary\n")
            fh.write("================================\n")
            fh.write(f"Identity threshold:       {IDENTITY_THRESHOLD*100:.0f}%\n")
            fh.write(f"Reference:                prokaryote train + val (rep_split1)\n")
            fh.write(f"Total eukaryote proteins: {total_proteins:,}\n")
            fh.write(f"Flagged (too similar):    {len(flagged_ids):,} ({len(flagged_ids)/total_proteins*100:.2f}%)\n")
            fh.write(f"Clean after filtering:    {total_proteins - len(flagged_ids):,}\n\n")
            fh.write("Breakdown by class:\n")
            for cls in sorted(total_by_class):
                n_total   = total_by_class[cls]
                n_flagged = flagged_by_class[cls]
                fh.write(f"  {cls:<20} {n_flagged:>6,} / {n_total:>6,}  ({n_flagged/n_total*100:.2f}%)\n")
        print(f"\nSummary saved: {OUT_SUMMARY}")

    print("\nDone.")


if __name__ == "__main__":
    main()
