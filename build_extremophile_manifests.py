from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

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
MANIFEST_BASE_FIELDS = [
    "assembly_accession",
    "taxon_name",
    "strain",
    "sequence_source_db",
    "sequence_source_accession",
    "sequence_level",
    "evidence_source",
    "local_fasta_file",
    "protein_count",
]
OPTIONAL_MANIFEST_FIELDS = [
    "psychrophile",
    "psychrophilic_alga",
    "thermophile",
    "hyperthermophile",
    "high_temperature_alga",
    "thermophilic_alga",
    "thermoacidophilic_alga",
    "ogt_celsius",
    "halophile",
    "saline_adapted_alga",
    "halophilic_alga",
    "halophile_class",
    "nacl_optimum_M",
    "acidophile",
    "acidophilic_alga",
    "acidophile_class",
    "alkaliphile",
    "alkaliphile_class",
    "ph_optimum",
]
MISSING_FIELDS = ["assembly_accession", "taxon_name", "strain"]

DATASETS = [
    {
        "label": "psychrophile",
        "tracker": BASE / "sequence_retrieval_tracker_psychrophile_yes.tsv",
        "fasta_dir": BASE / "psychrophile_raw_protein_fastas",
        "manifest_dir": BASE / "psychrophile_manifests",
        "merged_dir": BASE / "psychrophile_merged",
    },
    {
        "label": "thermophile",
        "tracker": BASE / "sequence_retrieval_tracker_tempura_bacdive_yes.tsv",
        "fasta_dir": BASE / "thermophile_raw_protein_fastas",
        "manifest_dir": BASE / "thermophile_manifests",
        "merged_dir": BASE / "thermophile_merged",
    },
    {
        "label": "halophile",
        "tracker": BASE / "sequence_retrieval_tracker_halophile_bacdive_ready.tsv",
        "fasta_dir": BASE / "halophile_raw_protein_fastas",
        "manifest_dir": BASE / "halophile_manifests",
        "merged_dir": BASE / "halophile_merged",
    },
    {
        "label": "acidophile",
        "tracker": BASE / "sequence_retrieval_tracker_acidophile_bacdive_yes.tsv",
        "fasta_dir": BASE / "acidophile_raw_protein_fastas",
        "manifest_dir": BASE / "acidophile_manifests",
        "merged_dir": BASE / "acidophile_merged",
    },
    {
        "label": "moderate_acidophile",
        "tracker": BASE / "sequence_retrieval_tracker_moderate_acidophile_bacdive_yes.tsv",
        "fasta_dir": BASE / "moderate_acidophile_raw_protein_fastas",
        "manifest_dir": BASE / "moderate_acidophile_manifests",
        "merged_dir": BASE / "moderate_acidophile_merged",
    },
    {
        "label": "acidotolerant",
        "tracker": BASE / "sequence_retrieval_tracker_acidotolerant_bacdive_yes.tsv",
        "fasta_dir": BASE / "acidotolerant_raw_protein_fastas",
        "manifest_dir": BASE / "acidotolerant_manifests",
        "merged_dir": BASE / "acidotolerant_merged",
    },
    {
        "label": "alkaliphile",
        "tracker": BASE / "sequence_retrieval_tracker_alkaliphile_bacdive_yes.tsv",
        "fasta_dir": BASE / "alkaliphile_raw_protein_fastas",
        "manifest_dir": BASE / "alkaliphile_manifests",
        "merged_dir": BASE / "alkaliphile_merged",
    },
    {
        "label": "eukaryotic_algae",
        "tracker": BASE / "sequence_retrieval_tracker_eukaryotic_algae_yes.tsv",
        "fasta_dir": BASE / "algae_raw_protein_fastas",
        "manifest_dir": BASE / "algae_manifests",
        "merged_dir": BASE / "algae_merged",
    },
]
COMBINED_DATASETS = [
    {
        "label": "acidophile_combined",
        "manifest_label": "combined_acidophiles",
        "members": [
            {
                "tracker": BASE / "sequence_retrieval_tracker_acidophile_bacdive_yes.tsv",
                "fasta_dir": BASE / "acidophile_raw_protein_fastas",
            },
            {
                "tracker": BASE / "sequence_retrieval_tracker_moderate_acidophile_bacdive_yes.tsv",
                "fasta_dir": BASE / "moderate_acidophile_raw_protein_fastas",
            },
        ],
        "manifest_dir": BASE / "acidophile_manifests",
        "merged_dir": BASE / "acidophile_merged",
    },
]


def build_manifest_names(
    label: str,
    matched: int,
    expected: int,
    manifest_label: Optional[str] = None,
) -> Tuple[str, str]:
    # When some FASTAs are still missing, the filename keeps both numbers so
    # you can see the gap immediately.
    manifest_stem = manifest_label or label
    if matched == expected:
        return (
            f"sequences_manifest_{manifest_stem}_{matched}.tsv",
            f"{label}_core_{matched}_merged.faa",
        )
    return (
        f"sequences_manifest_{manifest_stem}_{matched}_of_{expected}.tsv",
        f"{label}_core_{matched}_of_{expected}_merged.faa",
    )


def collect_manifest_rows(
    tracker_path: Path,
    fasta_dir: Path,
) -> Tuple[List[Dict[str, object]], List[Path], List[Dict[str, str]], List[str], int]:
    rows = read_tsv(tracker_path)
    optional_fields = [field for field in OPTIONAL_MANIFEST_FIELDS if rows and field in rows[0]]
    manifest_rows = []
    fasta_paths = []
    missing_rows = []

    for row in rows:
        accession = strip_accession_prefix(row["sequence_source_accession"])
        matches = find_matching_fasta_files(fasta_dir, accession)
        if not matches:
            missing_rows.append(
                {
                    "assembly_accession": accession,
                    "taxon_name": row["taxon_name"],
                    "strain": row["strain"],
                }
            )
            continue

        fasta_path = matches[0]
        fasta_paths.append(fasta_path)
        protein_count = count_fasta_records(fasta_path)

        manifest_row = {
            "assembly_accession": accession,
            "taxon_name": row["taxon_name"],
            "strain": row["strain"],
            "sequence_source_db": row["sequence_source_db"],
            "sequence_source_accession": row["sequence_source_accession"],
            "sequence_level": row["sequence_level"],
            "evidence_source": row["evidence_source"],
            "local_fasta_file": fasta_path.name,
            "protein_count": protein_count,
        }
        manifest_row.update({field: row.get(field, "") for field in optional_fields})

        manifest_rows.append(manifest_row)

    return manifest_rows, fasta_paths, missing_rows, optional_fields, len(rows)


def write_manifest_bundle(
    label: str,
    manifest_rows: List[Dict[str, object]],
    fasta_paths: List[Path],
    missing_rows: List[Dict[str, str]],
    optional_fields: List[str],
    expected: int,
    manifest_dir: Path,
    merged_dir: Path,
    manifest_label: Optional[str] = None,
) -> Dict[str, object]:
    manifest_dir.mkdir(exist_ok=True)
    merged_dir.mkdir(exist_ok=True)

    manifest_rows.sort(key=lambda item: (item["taxon_name"], item["strain"]))
    matched = len(manifest_rows)
    manifest_name, merged_name = build_manifest_names(label, matched, expected, manifest_label)

    manifest_path = manifest_dir / manifest_name
    merged_path = merged_dir / merged_name
    missing_path = manifest_dir / f"{label}_missing_fastas.tsv"
    manifest_fields = MANIFEST_BASE_FIELDS + optional_fields

    # The manifest is the structured summary; the merged FASTA is the actual
    # sequence bundle used downstream.
    write_tsv(manifest_path, manifest_rows, manifest_fields)
    merge_fasta_files(fasta_paths, merged_path)

    if missing_rows:
        write_tsv(missing_path, missing_rows, MISSING_FIELDS)
    elif missing_path.exists():
        missing_path.unlink()

    return {
        "label": label,
        "expected_rows": expected,
        "manifest_rows": matched,
        "manifest_path": manifest_path,
        "merged_path": merged_path,
        "missing_count": len(missing_rows),
        "missing_path": missing_path if missing_rows else None,
        "total_proteins": sum(row["protein_count"] for row in manifest_rows),
    }


def build_dataset(dataset: Dict[str, Any]) -> Dict[str, object]:
    # One pass over the tracker does three things:
    # 1. find the local FASTA,
    # 2. store metadata in a manifest row,
    # 3. remember the FASTA for the merged output file.
    label = str(dataset["label"])
    tracker_path = Path(dataset["tracker"])
    fasta_dir = Path(dataset["fasta_dir"])
    manifest_dir = Path(dataset["manifest_dir"])
    merged_dir = Path(dataset["merged_dir"])

    manifest_rows, fasta_paths, missing_rows, optional_fields, expected = collect_manifest_rows(
        tracker_path,
        fasta_dir,
    )
    return write_manifest_bundle(
        label,
        manifest_rows,
        fasta_paths,
        missing_rows,
        optional_fields,
        expected,
        manifest_dir,
        merged_dir,
        dataset.get("manifest_label"),
    )


def build_combined_dataset(dataset: Dict[str, Any]) -> Dict[str, object]:
    label = str(dataset["label"])
    manifest_dir = Path(dataset["manifest_dir"])
    merged_dir = Path(dataset["merged_dir"])

    manifest_rows: List[Dict[str, object]] = []
    fasta_paths: List[Path] = []
    missing_rows: List[Dict[str, str]] = []
    optional_fields: List[str] = []
    expected = 0

    for member in dataset["members"]:
        member_rows, member_fastas, member_missing, member_optional_fields, member_expected = collect_manifest_rows(
            Path(member["tracker"]),
            Path(member["fasta_dir"]),
        )
        manifest_rows.extend(member_rows)
        fasta_paths.extend(member_fastas)
        missing_rows.extend(member_missing)
        expected += member_expected
        for field in member_optional_fields:
            if field not in optional_fields:
                optional_fields.append(field)

    return write_manifest_bundle(
        label,
        manifest_rows,
        fasta_paths,
        missing_rows,
        optional_fields,
        expected,
        manifest_dir,
        merged_dir,
        dataset.get("manifest_label"),
    )


def main() -> None:
    for dataset in DATASETS:
        result = build_dataset(dataset)
        print(f"{result['label']}:")
        print(f"  manifest -> {result['manifest_path']}")
        print(f"  merged   -> {result['merged_path']}")
        print(f"  rows     -> {result['manifest_rows']} / {result['expected_rows']}")
        print(f"  proteins -> {result['total_proteins']}")
        print(f"  missing  -> {result['missing_count']}")
        if result["missing_path"]:
            print(f"  missing list -> {result['missing_path']}")
    for dataset in COMBINED_DATASETS:
        result = build_combined_dataset(dataset)
        print(f"{result['label']}:")
        print(f"  manifest -> {result['manifest_path']}")
        print(f"  merged   -> {result['merged_path']}")
        print(f"  rows     -> {result['manifest_rows']} / {result['expected_rows']}")
        print(f"  proteins -> {result['total_proteins']}")
        print(f"  missing  -> {result['missing_count']}")
        if result["missing_path"]:
            print(f"  missing list -> {result['missing_path']}")


if __name__ == "__main__":
    main()
