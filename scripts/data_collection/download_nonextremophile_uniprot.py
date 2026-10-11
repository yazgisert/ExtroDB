# Downloads all reviewed bacterial Swiss-Prot sequences from UniProt (taxonomy_id:2), removes the 247 organisms
# of the master index by NCBI taxon ID and writes a random sample of 142,000 (seed 42) as nonextremophile.faa.

import csv
import gzip
import io
import random
import time
import urllib.request
from pathlib import Path

BASE_DIR   = Path(__file__).resolve().parents[2]
INDEX_PATH = BASE_DIR / "searchable_resource" / "extremophile_master_index.tsv"
OUT_DIR    = BASE_DIR / "nonextremophile_data"
OUT_FAA    = OUT_DIR / "nonextremophile.faa"

# How many non-extremophile sequences to keep (matches our extremophile split size)
TARGET_N = 142_000

SEED = 42
random.seed(SEED)

UNIPROT_URL = (
    "https://rest.uniprot.org/uniprotkb/stream"
    "?query=(reviewed:true)%20AND%20(taxonomy_id:2)"
    "&format=fasta"
    "&compressed=true"
)


def load_excluded_taxids():
    with open(INDEX_PATH) as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    ids = set(r["ncbi_taxid"].strip() for r in rows if r.get("ncbi_taxid","").strip())
    print(f"Loaded {len(ids)} taxon IDs to exclude.")
    return ids


def parse_taxid(header):
    for part in header.split():
        if part.startswith("OX="):
            return part[3:].strip()
    return None


def parse_uniprot_fasta(raw_text, excluded_taxids):
    records = []
    current_id   = None
    current_desc = ""
    current_seq  = []
    current_skip = False

    for line in raw_text.splitlines():
        if line.startswith(">"):
            if current_id and not current_skip:
                records.append((current_id, current_desc, "".join(current_seq)))

            header  = line[1:]
            parts   = header.split("|")
            uniprot_id = parts[1] if len(parts) >= 2 else header.split()[0]
            desc       = header

            taxid = parse_taxid(header)
            current_skip = (taxid in excluded_taxids) if taxid else False

            current_id   = uniprot_id
            current_desc = desc
            current_seq  = []
        else:
            current_seq.append(line.strip())

    if current_id and not current_skip:
        records.append((current_id, current_desc, "".join(current_seq)))

    return records


def download_and_filter(excluded_taxids):
    print(f"Downloading Swiss-Prot reviewed bacterial sequences from UniProt ...")

    req = urllib.request.Request(UNIPROT_URL, headers={"User-Agent": "Python/thesis-download"})

    with urllib.request.urlopen(req, timeout=300) as response:
        compressed = response.read()

    print(f"Downloaded {len(compressed) / 1e6:.1f} MB compressed. Decompressing ...")
    raw_text = gzip.decompress(compressed).decode("utf-8")
    print(f"Decompressed: {len(raw_text) / 1e6:.1f} MB")

    print("Parsing and filtering ...")
    records = parse_uniprot_fasta(raw_text, excluded_taxids)
    print(f"Sequences after excluding our {len(excluded_taxids)} organisms: {len(records)}")

    return records


def write_faa(records, out_path, target_n):
    if len(records) > target_n:
        print(f"Sampling {target_n} sequences from {len(records)} ...")
        records = random.sample(records, target_n)
    else:
        print(f"Keeping all {len(records)} sequences (fewer than target {target_n}).")

    out_path.parent.mkdir(exist_ok=True)
    with open(out_path, "w") as f:
        for uniprot_id, desc, seq in records:
            # Our standard format: >proteome_id||protein_id description||label
            f.write(f">{uniprot_id}||{uniprot_id} {desc}||non_extremophile\n")
            for i in range(0, len(seq), 60):
                f.write(seq[i:i+60] + "\n")

    print(f"\nSaved {len(records)} sequences to {out_path}")


if __name__ == "__main__":
    excluded = load_excluded_taxids()
    records  = download_and_filter(excluded)

    if not records:
        print("ERROR: No sequences downloaded. Check your internet connection.")
    else:
        write_faa(records, OUT_FAA, TARGET_N)
