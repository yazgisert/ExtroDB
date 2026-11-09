# Removes the proteins that could not be embedded with ProtT5 (all longer than 6,500 residues; the biocentral API
# fails on them even one at a time) from the split FASTAs, so that all four representations are trained and tested
# on the same proteins. The original FASTAs are kept in dataset_splits/backup_before_exclusion/ and are always the
# input, so the script can be run again. The removed proteins are listed in dataset_splits/excluded_unembeddable_proteins.tsv.
# Run from the repository root (biocentral env): python scripts/dataset_splits/exclude_unembeddable_proteins.py

import shutil
from pathlib import Path

import h5py

SPLITS = Path("dataset_splits")
BACKUP = SPLITS / "backup_before_exclusion"
OUT_TSV = SPLITS / "excluded_unembeddable_proteins.tsv"
RUNS = {"split1_30pct": "rep_split1", "split2_30pct": "rep_split2",
        "split1_20pct": "rep_split1_20pct", "split2_20pct": "rep_split2_20pct"}


def read_records(path):
    records, header, seq = [], None, []
    for line in open(path):
        line = line.rstrip("\n")
        if line.startswith(">"):
            if header is not None:
                records.append((header, "".join(seq)))
            header, seq = line, []
        else:
            seq.append(line.strip())
    if header is not None:
        records.append((header, "".join(seq)))
    return records


def protein_id(header):
    return header[1:].split("||")[1].split()[0].strip()


if __name__ == "__main__":
    BACKUP.mkdir(exist_ok=True)
    removed = []
    for split, prefix in RUNS.items():
        with h5py.File(f"embeddings/{split}.h5", "r") as f:
            embedded = set(f.keys())
        for part in ("train", "val", "test"):
            name = f"{prefix}_{part}.faa"
            if not (BACKUP / name).exists():
                shutil.copy2(SPLITS / name, BACKUP / name)
            records = read_records(BACKUP / name)
            kept = [(h, s) for h, s in records if protein_id(h) in embedded]
            for h, s in records:
                if protein_id(h) not in embedded:
                    removed.append((protein_id(h), len(s), h.split("||")[-1], split, part))
            with open(SPLITS / name, "w") as f:
                for h, s in kept:
                    f.write(f"{h}\n{s}\n")
            print(f"{name}: {len(kept):,} of {len(records):,} kept")

    with open(OUT_TSV, "w") as f:
        f.write("protein_id\tlength\tlabel\tsplit\tset\n")
        for row in removed:
            f.write("\t".join(map(str, row)) + "\n")
    print(f"\n{len(removed)} entries ({len({r[0] for r in removed})} unique proteins) removed; listed in {OUT_TSV}")
