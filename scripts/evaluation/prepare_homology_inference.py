# Writes the leakage-filtered algae proteins (query) and all prokaryote extremophile proteins (reference) as FASTA and
# prints the MMseqs2 search command. The query file is also the algae input for the embedding evaluations.

import csv
from pathlib import Path

BASE_DIR      = Path(__file__).resolve().parents[2]
INDEX_PATH    = BASE_DIR / "searchable_resource" / "extremophile_master_index.tsv"
BLACKLIST     = BASE_DIR / "results" / "eukaryote_leakage" / "flagged_protein_ids.txt"
OUT_DIR       = BASE_DIR / "results" / "homology_inference"

QUERY_FAA     = OUT_DIR / "filtered_eukaryote_query.faa"
REF_FAA       = OUT_DIR / "prokaryote_extremophile_ref.faa"
HITS_OUT      = OUT_DIR / "homology_hits.tsv"
TMP_DIR       = OUT_DIR / "mmseqs_tmp"

OUT_DIR.mkdir(parents=True, exist_ok=True)


def load_blacklist():
    with open(BLACKLIST) as f:
        return {line.strip() for line in f if line.strip()}


def build_query_fasta(blacklist):
    rows = []
    with open(INDEX_PATH) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            if row.get("domain_scope") != "eukaryote":
                continue
            fasta_path = (BASE_DIR / row["local_fasta_path"].strip()) if row.get("local_fasta_path") else None
            if not fasta_path or not fasta_path.exists():
                print(f"  WARNING: missing FASTA for {row['taxon_name']} — skipping")
                continue
            rows.append(fasta_path)

    written = skipped = 0
    with open(QUERY_FAA, "w") as out:
        for fasta_path in rows:
            current_id = None
            buf = []
            with open(fasta_path) as fh:
                for line in fh:
                    if line.startswith(">"):
                        if current_id is not None and current_id not in blacklist:
                            out.writelines(buf)
                            written += 1
                        elif current_id is not None:
                            skipped += 1
                        current_id = line[1:].split()[0]
                        buf = [line]
                    else:
                        buf.append(line)
                if current_id is not None:
                    if current_id not in blacklist:
                        out.writelines(buf)
                        written += 1
                    else:
                        skipped += 1

    print(f"Query FASTA: {written:,} sequences written, {skipped:,} blacklisted sequences skipped")
    print(f"  -> {QUERY_FAA}")
    return written


def build_reference_fasta():
    written = 0
    with open(INDEX_PATH) as f, open(REF_FAA, "w") as out:
        for row in csv.DictReader(f, delimiter="\t"):
            if row.get("domain_scope") == "eukaryote":
                continue
            fasta_path = (BASE_DIR / row["local_fasta_path"].strip()) if row.get("local_fasta_path") else None
            if not fasta_path or not fasta_path.exists():
                print(f"  WARNING: missing FASTA for {row['taxon_name']} — skipping")
                continue
            with open(fasta_path) as fh:
                for line in fh:
                    out.write(line)
                    if line.startswith(">"):
                        written += 1

    print(f"Reference FASTA: {written:,} prokaryote extremophile sequences")
    print(f"  -> {REF_FAA}")
    return written


def print_mmseqs_command():
    print("\nRun this MMseqs2 command:")
    print(f"""
mmseqs easy-search \\
    {QUERY_FAA} \\
    {REF_FAA} \\
    {HITS_OUT} \\
    {TMP_DIR} \\
    --format-output "query,target,pident,evalue,bits" \\
    -s 7.5 \\
    --threads 8

Then run:
    python scripts/evaluation/eval_homology_inference.py
""")


if __name__ == "__main__":
    print("Loading blacklist ...")
    blacklist = load_blacklist()
    print(f"  {len(blacklist):,} proteins blacklisted (>=30% identity to prokaryote train+val)\n")

    print("Building filtered eukaryote query FASTA ...")
    build_query_fasta(blacklist)

    print("\nBuilding prokaryote extremophile reference FASTA ...")
    build_reference_fasta()

    print_mmseqs_command()
