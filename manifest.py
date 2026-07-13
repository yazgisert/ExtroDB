from __future__ import annotations

from pathlib import Path

from common import (
    BASE_DIR,
    count_fasta_records,
    find_matching_fasta_files,
    merge_fasta_files,
    read_tsv,
    strip_accession_prefix,
    write_tsv,
)


BASE = BASE_DIR
TRACKER_FILE = BASE / "sequence_retrieval_tracker_tempura_bacdive_yes.tsv"
FASTA_DIR = BASE / "raw_protein_fastas"
MANIFEST_FILE = BASE / "sequences_manifest_tempura_33.tsv"
MERGED_FASTA = BASE / "thermophile_core_33_merged.faa"
FIELDNAMES = [
    "assembly_accession",
    "taxon_name",
    "strain",
    "thermophile",
    "hyperthermophile",
    "ogt_celsius",
    "local_fasta_file",
    "protein_count",
]


def build_manifest_rows(tracker_rows: list[dict[str, str]]) -> tuple[list[dict[str, object]], list[Path]]:
    # This is the older single-dataset manifest builder kept for the original
    # thermophile workflow.
    manifest_rows: list[dict[str, object]] = []
    fasta_paths: list[Path] = []

    for row in tracker_rows:
        accession = strip_accession_prefix(row["sequence_source_accession"])
        matches = find_matching_fasta_files(FASTA_DIR, accession)
        if not matches:
            print(f"Missing FASTA for {accession}")
            continue

        fasta_path = matches[0]
        fasta_paths.append(fasta_path)
        manifest_rows.append(
            {
                "assembly_accession": accession,
                "taxon_name": row["taxon_name"],
                "strain": row["strain"],
                "thermophile": row["thermophile"],
                "hyperthermophile": row["hyperthermophile"],
                "ogt_celsius": row["ogt_celsius"],
                "local_fasta_file": fasta_path.name,
                "protein_count": count_fasta_records(fasta_path),
            }
        )

    return manifest_rows, fasta_paths


def main() -> None:
    tracker_rows = read_tsv(TRACKER_FILE)
    manifest_rows, fasta_paths = build_manifest_rows(tracker_rows)
    write_tsv(MANIFEST_FILE, manifest_rows, FIELDNAMES)
    merge_fasta_files(fasta_paths, MERGED_FASTA)

    print("Done")
    print("Manifest:", MANIFEST_FILE)
    print("Merged FASTA:", MERGED_FASTA)
    print("Rows:", len(manifest_rows))


if __name__ == "__main__":
    main()
