# Combines the protein FASTAs of all single-label prokaryotes into one file for MMseqs2 clustering (headers become
# >proteome_id||original_header) and writes one metadata row per organism with the 5-class label (phenotype_focus,
# already collapsed in build_searchable_extremophile_index.py) and the 12-class label.

from __future__ import annotations

import csv
from pathlib import Path


BASE_DIR   = Path(__file__).resolve().parents[2]
SPLITS_DIR = BASE_DIR / "dataset_splits"

MASTER_INDEX = BASE_DIR / "searchable_resource" / "extremophile_master_index.tsv"
OUT_FASTA    = SPLITS_DIR / "split_input_prokaryote_singlelabel.faa"
OUT_META     = SPLITS_DIR / "split_input_metadata.tsv"

VALID_SPLIT1_LABELS = {"thermophile", "halophile", "acidophile", "alkaliphile", "psychrophile"}


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def is_single_label_prokaryote(row: dict[str, str]) -> bool:
    return (
        row.get("domain_scope") == "prokaryote"
        and "|" not in row.get("phenotype_focus", "")
        and row.get("proteome_status") == "downloaded"
        and row.get("local_fasta_path", "").strip() != ""
    )


def derive_split2_label(row: dict[str, str]) -> str:
    # 12-class label from phenotype_labels; within each class the most specific subtype is checked first.
    labels_raw = row.get("phenotype_labels", "").lower()

    if "acidotolerant" in labels_raw:
        return "acidotolerant"
    if "moderate acidophile" in labels_raw:
        return "moderate_acidophile"
    if "acidophile" in labels_raw:
        return "acidophile"

    # Halophile subclasses — format is "halophile (slight)" etc.
    if "halophile (slight)" in labels_raw:
        return "slight_halophile"
    if "halophile (moderate)" in labels_raw:
        return "moderate_halophile"
    if "halophile (extreme)" in labels_raw:
        return "extreme_halophile"
    if "halophile" in labels_raw:
        return "halophile"          # fallback for any future plain halophile entry

    if "alkalitolerant" in labels_raw:
        return "alkalitolerant"
    if "extreme alkaliphile" in labels_raw:
        return "extreme_alkaliphile"
    if "alkaliphile" in labels_raw:
        return "alkaliphile"

    if "hyperthermophile" in labels_raw:
        return "hyperthermophile"
    if "thermophile" in labels_raw:
        return "thermophile"

    if "psychrophile" in labels_raw:
        return "psychrophile"

    # Should never reach here for a valid single-label prokaryote
    raise ValueError(
        f"Could not derive Split 2 label for {row.get('proteome_id')} "
        f"(phenotype_labels: {row.get('phenotype_labels')!r})"
    )


def main() -> None:
    print(f"Reading master index: {MASTER_INDEX}")
    all_rows = read_tsv(MASTER_INDEX)

    organisms = [r for r in all_rows if is_single_label_prokaryote(r)]
    print(f"Single-label prokaryotes (downloaded): {len(organisms)}")

    # Safety assertion: confirm every organism has a valid Split 1 label
    bad_labels = [r for r in organisms if r.get("phenotype_focus") not in VALID_SPLIT1_LABELS]
    if bad_labels:
        raise ValueError(
            f"{len(bad_labels)} organism(s) have unexpected phenotype_focus values:\n"
            + "\n".join(f"  {r['proteome_id']}: {r['phenotype_focus']}" for r in bad_labels)
        )

    missing = []
    for row in organisms:
        fasta_path = BASE_DIR / row["local_fasta_path"]
        if not fasta_path.exists():
            missing.append((row["proteome_id"], str(fasta_path)))

    if missing:
        print(f"\nWARNING: {len(missing)} FASTA file(s) not found:")
        for pid, path in missing:
            print(f"  {pid}: {path}")
        print("These organisms will be skipped.\n")
        organisms = [r for r in organisms
                     if (BASE_DIR / r["local_fasta_path"]).exists()]

    print(f"\nWriting combined FASTA: {OUT_FASTA}")
    total_proteins = 0
    with OUT_FASTA.open("w", encoding="utf-8") as out_fh:
        for row in organisms:
            proteome_id = row["proteome_id"]
            fasta_path  = BASE_DIR / row["local_fasta_path"]
            with fasta_path.open(encoding="utf-8", errors="replace") as fasta_fh:
                for line in fasta_fh:
                    line = line.rstrip("\n")
                    if line.startswith(">"):
                        original = line[1:]
                        out_fh.write(f">{proteome_id}||{original}\n")
                        total_proteins += 1
                    else:
                        out_fh.write(line + "\n")

    print(f"  Total proteins written: {total_proteins:,}")

    META_FIELDS = [
        "proteome_id",
        "taxon_name",
        "strain",
        "phenotype_focus",    # Split 1 label (coarse, 5 classes)
        "split2_label",       # Split 2 label (subtype)
        "phenotype_labels",   # Full label string from master index
        "ogt_celsius",
        "nacl_optimum_M",
        "ph_optimum",
        "sequence_source_accession",
        "protein_count",
        "local_fasta_path",
    ]

    print(f"Writing metadata TSV:   {OUT_META}")
    with OUT_META.open("w", encoding="utf-8", newline="") as meta_fh:
        writer = csv.DictWriter(meta_fh, fieldnames=META_FIELDS, delimiter="\t",
                                extrasaction="ignore")
        writer.writeheader()
        for row in organisms:
            out_row = {f: row.get(f, "") for f in META_FIELDS}
            out_row["split2_label"] = derive_split2_label(row)
            writer.writerow(out_row)

    print(f"  Organisms written: {len(organisms)}")
    print("\nDone. Next step: run MMseqs2 easy-cluster on the combined FASTA.")
    print(f"  Input:  {OUT_FASTA}")


if __name__ == "__main__":
    main()
