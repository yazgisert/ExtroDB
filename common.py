from __future__ import annotations

import csv
import re
from pathlib import Path
from typing import Iterable, Mapping


BASE_DIR = Path(__file__).resolve().parent
TSV_DELIMITER = "\t"


def read_tsv(path: Path) -> list[dict[str, str]]:
    # Most scripts in this folder read small TSV tables into memory, so a
    # shared helper keeps the calling code short and consistent.
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter=TSV_DELIMITER))


def write_tsv(
    path: Path,
    rows: Iterable[Mapping[str, object]],
    fieldnames: Iterable[str] | None = None,
) -> None:
    row_list = [dict(row) for row in rows]
    columns = list(fieldnames or (row_list[0].keys() if row_list else ()))
    if not columns:
        raise ValueError("fieldnames are required when writing an empty TSV.")

    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter=TSV_DELIMITER)
        writer.writeheader()
        writer.writerows(row_list)


def strip_accession_prefix(value: str) -> str:
    # GTDB-style accessions sometimes carry RS_/GB_ prefixes that do not match
    # local FASTA filenames, so we normalize them once here.
    return value[3:] if value.startswith(("RS_", "GB_")) else value


def parse_ncbi_assembly_accession(value: str) -> tuple[str, str, str] | None:
    normalized = strip_accession_prefix(value.strip())
    match = re.fullmatch(r"(GC[AF])_(\d+)(?:\.(\d+))?", normalized)
    if not match:
        return None
    prefix, number, version = match.groups()
    return prefix, number, version or ""


def accessions_equivalent(left: str, right: str) -> bool:
    left_normalized = strip_accession_prefix(left.strip())
    right_normalized = strip_accession_prefix(right.strip())
    if left_normalized == right_normalized:
        return True

    left_parts = parse_ncbi_assembly_accession(left_normalized)
    right_parts = parse_ncbi_assembly_accession(right_normalized)
    if left_parts is None or right_parts is None:
        return False

    _, left_number, left_version = left_parts
    _, right_number, right_version = right_parts
    if left_number != right_number:
        return False
    if left_version and right_version and left_version != right_version:
        return False
    return True


def find_matching_fasta_files(fasta_dir: Path, accession: str) -> list[Path]:
    # Accept either GenBank/RefSeq assembly counterparts and tolerate missing
    # version suffixes in the tracker because NCBI downloads usually include
    # the explicit version in the filename.
    matches = []
    for fasta_path in sorted(fasta_dir.glob("*.faa")):
        local_accession = fasta_path.name.split("__", 1)[0]
        if accessions_equivalent(accession, local_accession):
            matches.append(fasta_path)
    return matches


def build_taxon_key(row: Mapping[str, str]) -> tuple[str, str]:
    # In this project, taxon name + strain is the main join key between
    # curation tables, trackers, and manifests.
    return (row.get("taxon_name", "").strip(), row.get("strain", "").strip())


def count_fasta_records(path: Path) -> int:
    # Each FASTA record starts with ">", so counting those lines gives the
    # protein count stored in manifests and summaries.
    with path.open("r", encoding="utf-8", errors="ignore") as handle:
        return sum(1 for line in handle if line.startswith(">"))


def merge_fasta_files(paths: Iterable[Path], destination: Path) -> None:
    # The merged FASTA is just the ordered concatenation of all selected
    # proteomes for a dataset branch.
    with destination.open("w", encoding="utf-8") as out_handle:
        for fasta_path in sorted(paths):
            with fasta_path.open("r", encoding="utf-8", errors="ignore") as in_handle:
                out_handle.write(in_handle.read())
            out_handle.write("\n")


def relative_to_base(path: Path | None, base: Path = BASE_DIR) -> str:
    # Relative paths make the output tables portable and easier to read.
    if path is None:
        return ""
    try:
        return str(path.relative_to(base))
    except ValueError:
        return str(path)
