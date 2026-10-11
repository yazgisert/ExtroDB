# Builds the master index of all curated proteomes from the curation tables, retrieval trackers and manifests of
# every class: one membership row per (proteome, class) and one master row per proteome with all labels, measured
# values (OGT, NaCl, pH) and the local FASTA path. Writes the master index, the membership index and a summary to
# searchable_resource/.

from __future__ import annotations

import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from common import (
    BASE_DIR,
    build_taxon_key,
    read_tsv,
    relative_to_base,
    strip_accession_prefix,
    write_tsv,
)


BASE = BASE_DIR
OUTPUT_DIR = BASE / "searchable_resource"
MULTI_VALUE_SEPARATOR = "|"
TEXT_LIST_SEPARATOR = " || "

RAW_MEMBERSHIP_FIELDNAMES = [
    "dataset_branch",
    "dataset_layer",
    "domain_scope",
    "phenotype_focus",
    "taxon_name",
    "ncbi_taxid",
    "strain",
    "strain_aliases",
    "identity_level",
    "collection_source",
    "collection_match",
    "psychrophile",
    "psychrophilic_alga",
    "thermophile",
    "hyperthermophile",
    "high_temperature_alga",
    "thermophilic_alga",
    "thermoacidophilic_alga",
    "halophile",
    "saline_adapted_alga",
    "halophilic_alga",
    "halophile_class",
    "acidophile",
    "acidophilic_alga",
    "acidophile_class",
    "alkaliphile",
    "alkaliphile_class",
    "ogt_celsius",
    "nacl_optimum_M",
    "ph_optimum",
    "evidence_type",
    "evidence_source",
    "include_in_core",
    "proteome_status",
    "proteome_in_manifest",
    "protein_count",
    "sequence_source_db",
    "sequence_source_accession",
    "sequence_level",
    "retrieval_status",
    "retrieval_notes",
    "source_table_file",
    "retrieval_tracker_file",
    "manifest_file",
    "merged_fasta_file",
    "local_fasta_file",
    "local_fasta_path",
    "notes",
]
MEMBERSHIP_FIELDNAMES = ["proteome_id"] + RAW_MEMBERSHIP_FIELDNAMES
MASTER_FIELDNAMES = [
    "proteome_id",
    "dataset_branch",
    "dataset_layer",
    "domain_scope",
    "phenotype_focus",
    "phenotype_labels",
    "branch_membership_count",
    "taxon_name",
    "ncbi_taxid",
    "strain",
    "strain_aliases",
    "identity_level",
    "collection_source",
    "collection_match",
    "psychrophile",
    "psychrophilic_alga",
    "thermophile",
    "hyperthermophile",
    "high_temperature_alga",
    "thermophilic_alga",
    "thermoacidophilic_alga",
    "halophile",
    "saline_adapted_alga",
    "halophilic_alga",
    "halophile_class",
    "acidophile",
    "acidophilic_alga",
    "acidophile_class",
    "alkaliphile",
    "alkaliphile_class",
    "ogt_celsius",
    "nacl_optimum_M",
    "ph_optimum",
    "evidence_type",
    "evidence_source",
    "include_in_core",
    "proteome_status",
    "proteome_in_manifest",
    "protein_count",
    "sequence_source_db",
    "sequence_source_accession",
    "sequence_level",
    "retrieval_status",
    "retrieval_notes",
    "source_table_file",
    "retrieval_tracker_file",
    "manifest_file",
    "merged_fasta_file",
    "local_fasta_file",
    "local_fasta_path",
    "local_fasta_paths",
    "notes",
]
SUMMARY_FIELDNAMES = [
    "dataset_branch",
    "dataset_layer",
    "proteome_status",
    "row_count",
    "protein_count_sum",
]

DATASETS = [
    {
        "dataset_branch": "psychrophile_core",
        "dataset_layer": "core",
        "domain_scope": "prokaryote",
        "phenotype_focus": "psychrophile",
        "source_table": BASE / "curation" / "core_taxon_tables" / "core_taxon_table_psychrophile.tsv",
        "retrieval_tracker": BASE / "curation" / "retrieval_trackers" / "sequence_retrieval_tracker_psychrophile_yes.tsv",
        "manifest_tracker": BASE / "curation" / "retrieval_trackers" / "sequence_retrieval_tracker_psychrophile_yes.tsv",
        "manifest_dir": BASE / "curation" / "manifests" / "psychrophile_manifests",
        "merged_dir": BASE / "psychrophile_merged",
        "raw_fasta_dir": BASE / "psychrophile_raw_protein_fastas",
    },
    {
        "dataset_branch": "thermophile_core",
        "dataset_layer": "core",
        "domain_scope": "prokaryote",
        "phenotype_focus": "thermophile",
        "source_table": BASE / "curation" / "core_taxon_tables" / "core_taxon_table_tempura_bacdive_yes.tsv",
        "retrieval_tracker": BASE / "curation" / "retrieval_trackers" / "sequence_retrieval_tracker_tempura_bacdive_yes.tsv",
        "manifest_tracker": BASE / "curation" / "retrieval_trackers" / "sequence_retrieval_tracker_tempura_bacdive_yes.tsv",
        "manifest_dir": BASE / "curation" / "manifests" / "thermophile_manifests",
        "merged_dir": BASE / "thermophile_merged",
        "raw_fasta_dir": BASE / "thermophile_raw_protein_fastas",
    },
    {
        "dataset_branch": "halophile_core",
        "dataset_layer": "core",
        "domain_scope": "prokaryote",
        "phenotype_focus": "halophile",
        "source_table": BASE / "curation" / "core_taxon_tables" / "core_taxon_table_halophile_bacdive.tsv",
        "retrieval_tracker": BASE / "curation" / "retrieval_trackers" / "sequence_retrieval_tracker_halophile_bacdive_yes.tsv",
        "manifest_tracker": BASE / "curation" / "retrieval_trackers" / "sequence_retrieval_tracker_halophile_bacdive_ready.tsv",
        "manifest_dir": BASE / "curation" / "manifests" / "halophile_manifests",
        "merged_dir": BASE / "halophile_merged",
        "raw_fasta_dir": BASE / "halophile_raw_protein_fastas",
    },
    {
        "dataset_branch": "eukaryotic_algae_exploratory",
        "dataset_layer": "exploratory",
        "domain_scope": "eukaryote",
        "phenotype_focus": "eukaryotic_algae",
        "source_table": BASE / "curation" / "core_taxon_tables" / "core_taxon_table_eukaryotic_algae.tsv",
        "retrieval_tracker": BASE / "curation" / "retrieval_trackers" / "sequence_retrieval_tracker_eukaryotic_algae_yes.tsv",
        "manifest_tracker": BASE / "curation" / "retrieval_trackers" / "sequence_retrieval_tracker_eukaryotic_algae_yes.tsv",
        "manifest_dir": BASE / "curation" / "manifests" / "algae_manifests",
        "merged_dir": BASE / "algae_merged",
        "raw_fasta_dir": BASE / "algae_raw_protein_fastas",
    },
    {
        "dataset_branch": "acidophile_expansion",
        "dataset_layer": "expansion",
        "domain_scope": "prokaryote",
        "phenotype_focus": "acidophile",
        "source_table": BASE / "curation" / "core_taxon_tables" / "core_taxon_table_acidophile_bacdive.tsv",
        "retrieval_tracker": BASE / "curation" / "retrieval_trackers" / "sequence_retrieval_tracker_acidophile_bacdive_yes.tsv",
        "manifest_tracker": BASE / "curation" / "retrieval_trackers" / "sequence_retrieval_tracker_acidophile_bacdive_yes.tsv",
        "manifest_dir": BASE / "curation" / "manifests" / "acidophile_manifests",
        "merged_dir": BASE / "acidophile_merged",
        "raw_fasta_dir": BASE / "acidophile_raw_protein_fastas",
    },
    {
        "dataset_branch": "moderate_acidophile_expansion",
        "dataset_layer": "expansion",
        "domain_scope": "prokaryote",
        "phenotype_focus": "moderate_acidophile",
        "source_table": BASE / "curation" / "core_taxon_tables" / "core_taxon_table_moderate_acidophile_bacdive.tsv",
        "retrieval_tracker": BASE / "curation" / "retrieval_trackers" / "sequence_retrieval_tracker_moderate_acidophile_bacdive_yes.tsv",
        "manifest_tracker": BASE / "curation" / "retrieval_trackers" / "sequence_retrieval_tracker_moderate_acidophile_bacdive_yes.tsv",
        "manifest_dir": BASE / "curation" / "manifests" / "moderate_acidophile_manifests",
        "merged_dir": BASE / "moderate_acidophile_merged",
        "raw_fasta_dir": BASE / "moderate_acidophile_raw_protein_fastas",
    },
    {
        "dataset_branch": "acidotolerant_expansion",
        "dataset_layer": "expansion",
        "domain_scope": "prokaryote",
        "phenotype_focus": "acidotolerant",
        "source_table": BASE / "curation" / "core_taxon_tables" / "core_taxon_table_acidotolerant_bacdive.tsv",
        "retrieval_tracker": BASE / "curation" / "retrieval_trackers" / "sequence_retrieval_tracker_acidotolerant_bacdive_yes.tsv",
        "manifest_tracker": BASE / "curation" / "retrieval_trackers" / "sequence_retrieval_tracker_acidotolerant_bacdive_yes.tsv",
        "manifest_dir": BASE / "curation" / "manifests" / "acidotolerant_manifests",
        "merged_dir": BASE / "acidotolerant_merged",
        "raw_fasta_dir": BASE / "acidotolerant_raw_protein_fastas",
    },
    {
        "dataset_branch": "alkaliphile_core",
        "dataset_layer": "core",
        "domain_scope": "prokaryote",
        "phenotype_focus": "alkaliphile",
        "source_table": BASE / "curation" / "core_taxon_tables" / "core_taxon_table_alkaliphile_bacdive.tsv",
        "retrieval_tracker": BASE / "curation" / "retrieval_trackers" / "sequence_retrieval_tracker_alkaliphile_bacdive_yes.tsv",
        "manifest_tracker": BASE / "curation" / "retrieval_trackers" / "sequence_retrieval_tracker_alkaliphile_bacdive_yes.tsv",
        "manifest_dir": BASE / "curation" / "manifests" / "alkaliphile_manifests",
        "merged_dir": BASE / "alkaliphile_merged",
        "raw_fasta_dir": BASE / "alkaliphile_raw_protein_fastas",
    },
]

MASTER_INDEX_PATH = OUTPUT_DIR / "extremophile_master_index.tsv"
MEMBERSHIP_INDEX_PATH = OUTPUT_DIR / "extremophile_membership_index.tsv"
SUMMARY_PATH = OUTPUT_DIR / "extremophile_master_summary.tsv"

STATUS_PRIORITY = {
    "downloaded": 5,
    "ready_no_local_fasta": 4,
    "ready": 3,
    "missing": 2,
    "tracker_present_status_blank": 1,
    "not_started": 0,
}
DATASET_LAYER_PRIORITY = {"core": 0, "expansion": 1, "exploratory": 2}


def resolve_manifest_and_merged_paths(dataset: Dict[str, object]) -> Tuple[Optional[Path], Optional[Path]]:
    tracker_path = dataset.get("manifest_tracker")
    manifest_dir = dataset.get("manifest_dir")
    merged_dir = dataset.get("merged_dir")
    if not isinstance(tracker_path, Path) or not isinstance(manifest_dir, Path) or not isinstance(merged_dir, Path):
        return None, None

    rows = read_tsv(tracker_path)
    expected = len(rows)
    label = str(dataset["phenotype_focus"])
    manifest_path = manifest_dir / f"sequences_manifest_{label}_{expected}.tsv"
    merged_path = merged_dir / f"{label}_core_{expected}_merged.faa"
    if manifest_path.exists() and merged_path.exists():
        return manifest_path, merged_path

    manifest_candidates = sorted(manifest_dir.glob(f"sequences_manifest_{label}_*_of_{expected}.tsv"))
    merged_candidates = sorted(merged_dir.glob(f"{label}_core_*_of_{expected}_merged.faa"))
    if manifest_candidates and merged_candidates:
        return manifest_candidates[-1], merged_candidates[-1]

    return manifest_path, merged_path


def build_lookup(rows: Iterable[Dict[str, str]]) -> Dict[Tuple[str, str], Dict[str, str]]:
    return {build_taxon_key(row): row for row in rows}


def normalize_join_token(value: str) -> str:
    return "".join(character.lower() for character in value if character.isalnum())


def relaxed_taxon_key(row: Dict[str, str]) -> Tuple[str, str]:
    return (
        normalize_join_token(row.get("taxon_name", "")),
        normalize_join_token(row.get("strain", "")),
    )


def build_relaxed_lookup(rows: Iterable[Dict[str, str]]) -> Dict[Tuple[str, str], Dict[str, str]]:
    lookup: Dict[Tuple[str, str], Dict[str, str]] = {}
    for row in rows:
        lookup.setdefault(relaxed_taxon_key(row), row)
    return lookup


def normalized_taxon_name(row: Dict[str, str]) -> str:
    return normalize_join_token(row.get("taxon_name", ""))


def build_unique_taxon_lookup(rows: Iterable[Dict[str, str]]) -> Dict[str, Dict[str, str]]:
    grouped: Dict[str, List[Dict[str, str]]] = defaultdict(list)
    for row in rows:
        grouped[normalized_taxon_name(row)].append(row)
    return {taxon: members[0] for taxon, members in grouped.items() if len(members) == 1}


def row_value(
    primary: Optional[Dict[str, str]],
    secondary: Optional[Dict[str, str]],
    field: str,
) -> str:
    if primary:
        value = str(primary.get(field, "") or "").strip()
        if value:
            return value
    if secondary:
        return str(secondary.get(field, "") or "").strip()
    return ""


def determine_proteome_status(
    tracker_row: Optional[Dict[str, str]],
    manifest_row: Optional[Dict[str, str]],
) -> str:
    if manifest_row:
        return "downloaded"
    if tracker_row is None:
        return "not_started"

    tracker_status = tracker_row.get("retrieval_status", "").strip()
    if tracker_status == "ready":
        return "ready_no_local_fasta"
    if tracker_status:
        return tracker_status
    return "tracker_present_status_blank"


def normalize_multi_values(values: Iterable[str]) -> List[str]:
    seen = []
    for value in values:
        cleaned = value.strip()
        if cleaned and cleaned not in seen:
            seen.append(cleaned)
    return seen


def join_multi_values(values: Iterable[str]) -> str:
    return MULTI_VALUE_SEPARATOR.join(normalize_multi_values(values))


def join_text_blobs(values: Iterable[str]) -> str:
    return TEXT_LIST_SEPARATOR.join(normalize_multi_values(values))


def split_multi_value(value: str) -> List[str]:
    return [item.strip() for item in value.split(MULTI_VALUE_SEPARATOR) if item.strip()]


def title_case_label(value: str) -> str:
    return " ".join(part.capitalize() for part in value.split("_"))


def choose_best_status(rows: List[Dict[str, object]]) -> str:
    return max(
        (str(row.get("proteome_status", "")) for row in rows),
        key=lambda status: (STATUS_PRIORITY.get(status, -1), status),
        default="",
    )


def sort_membership_rows(rows: List[Dict[str, object]]) -> List[Dict[str, object]]:
    return sorted(
        rows,
        key=lambda row: (
            -STATUS_PRIORITY.get(str(row.get("proteome_status", "")), -1),
            DATASET_LAYER_PRIORITY.get(str(row.get("dataset_layer", "")), 99),
            str(row.get("dataset_branch", "")),
            str(row.get("taxon_name", "")),
            str(row.get("strain", "")),
        ),
    )


def first_non_empty(rows: List[Dict[str, object]], field: str) -> str:
    for row in rows:
        value = str(row.get(field, "")).strip()
        if value:
            return value
    return ""


def unique_non_empty(rows: List[Dict[str, object]], field: str) -> List[str]:
    return normalize_multi_values(str(row.get(field, "")).strip() for row in rows)


def choose_sequence_accession(rows: List[Dict[str, object]]) -> str:
    values = unique_non_empty(rows, "sequence_source_accession")
    if not values:
        return ""
    preferred = sorted(values, key=lambda value: (strip_accession_prefix(value), value.startswith("GCF_"), value), reverse=True)
    return preferred[0]


def choose_canonical_local_path(rows: List[Dict[str, object]]) -> Tuple[str, str]:
    paths = unique_non_empty(rows, "local_fasta_path")
    files = unique_non_empty(rows, "local_fasta_file")
    if not paths:
        return (files[0] if files else "", "")
    canonical_path = sorted(paths)[0]
    canonical_file = Path(canonical_path).name
    return canonical_file, canonical_path


# Proteome ID = assembly accession, else FASTA file name, else taxon__strain.
def compute_proteome_id(row: Dict[str, object]) -> str:
    accession = strip_accession_prefix(str(row.get("sequence_source_accession", "")).strip())
    if accession:
        return accession
    local_fasta_path = str(row.get("local_fasta_path", "")).strip()
    if local_fasta_path:
        return Path(local_fasta_path).name
    taxon = str(row.get("taxon_name", "")).strip().replace(" ", "_")
    strain = str(row.get("strain", "")).strip().replace(" ", "_")
    return f"{taxon}__{strain}"


def build_membership_row(
    dataset: Dict[str, object],
    curation_row: Optional[Dict[str, str]],
    tracker_row: Optional[Dict[str, str]],
    manifest_row: Optional[Dict[str, str]],
) -> Dict[str, object]:
    primary_row = tracker_row or curation_row or {}
    secondary_row = curation_row if tracker_row else None
    local_fasta_file = manifest_row.get("local_fasta_file", "") if manifest_row else ""
    raw_fasta_dir = dataset.get("raw_fasta_dir")
    raw_fasta_path = ""
    if local_fasta_file and isinstance(raw_fasta_dir, Path):
        raw_fasta_path = relative_to_base(raw_fasta_dir / local_fasta_file)

    protein_count = manifest_row.get("protein_count", "") if manifest_row else ""
    sequence_source_db = tracker_row.get("sequence_source_db", "") if tracker_row else ""
    sequence_source_accession = tracker_row.get("sequence_source_accession", "") if tracker_row else ""
    sequence_level = tracker_row.get("sequence_level", "") if tracker_row else ""
    if manifest_row:
        sequence_source_db = sequence_source_db or manifest_row.get("sequence_source_db", "")
        sequence_source_accession = sequence_source_accession or manifest_row.get(
            "sequence_source_accession", ""
        )
        sequence_level = sequence_level or manifest_row.get("sequence_level", "")

    row = {
        "dataset_branch": dataset["dataset_branch"],
        "dataset_layer": dataset["dataset_layer"],
        "domain_scope": dataset["domain_scope"],
        "phenotype_focus": dataset["phenotype_focus"],
        "taxon_name": row_value(primary_row, secondary_row, "taxon_name"),
        "ncbi_taxid": row_value(primary_row, secondary_row, "ncbi_taxid"),
        "strain": row_value(primary_row, secondary_row, "strain"),
        "strain_aliases": row_value(primary_row, secondary_row, "strain_aliases"),
        "identity_level": row_value(primary_row, secondary_row, "identity_level"),
        "collection_source": row_value(primary_row, secondary_row, "collection_source"),
        "collection_match": row_value(primary_row, secondary_row, "collection_match"),
        "psychrophile": row_value(primary_row, secondary_row, "psychrophile"),
        "psychrophilic_alga": row_value(primary_row, secondary_row, "psychrophilic_alga"),
        "thermophile": row_value(primary_row, secondary_row, "thermophile"),
        "hyperthermophile": row_value(primary_row, secondary_row, "hyperthermophile"),
        "high_temperature_alga": row_value(primary_row, secondary_row, "high_temperature_alga"),
        "thermophilic_alga": row_value(primary_row, secondary_row, "thermophilic_alga"),
        "thermoacidophilic_alga": row_value(primary_row, secondary_row, "thermoacidophilic_alga"),
        "halophile": row_value(primary_row, secondary_row, "halophile"),
        "saline_adapted_alga": row_value(primary_row, secondary_row, "saline_adapted_alga"),
        "halophilic_alga": row_value(primary_row, secondary_row, "halophilic_alga"),
        "halophile_class": row_value(primary_row, secondary_row, "halophile_class"),
        "acidophile": row_value(primary_row, secondary_row, "acidophile"),
        "acidophilic_alga": row_value(primary_row, secondary_row, "acidophilic_alga"),
        "acidophile_class": row_value(primary_row, secondary_row, "acidophile_class"),
        "alkaliphile": row_value(primary_row, secondary_row, "alkaliphile"),
        "alkaliphile_class": row_value(primary_row, secondary_row, "alkaliphile_class"),
        "ogt_celsius": row_value(primary_row, secondary_row, "ogt_celsius"),
        "nacl_optimum_M": row_value(primary_row, secondary_row, "nacl_optimum_M"),
        "ph_optimum": row_value(primary_row, secondary_row, "ph_optimum"),
        "evidence_type": row_value(primary_row, secondary_row, "evidence_type"),
        "evidence_source": row_value(primary_row, secondary_row, "evidence_source"),
        "include_in_core": row_value(primary_row, secondary_row, "include_in_core"),
        "proteome_status": determine_proteome_status(tracker_row, manifest_row),
        "proteome_in_manifest": "yes" if manifest_row else "no",
        "protein_count": protein_count,
        "sequence_source_db": sequence_source_db,
        "sequence_source_accession": sequence_source_accession,
        "sequence_level": sequence_level,
        "retrieval_status": tracker_row.get("retrieval_status", "") if tracker_row else "",
        "retrieval_notes": tracker_row.get("retrieval_notes", "") if tracker_row else "",
        "source_table_file": relative_to_base(dataset["source_table"]),
        "retrieval_tracker_file": relative_to_base(dataset.get("retrieval_tracker")),
        "manifest_file": relative_to_base(dataset.get("manifest")),
        "merged_fasta_file": relative_to_base(dataset.get("merged_fasta")),
        "local_fasta_file": local_fasta_file,
        "local_fasta_path": raw_fasta_path,
        "notes": row_value(primary_row, secondary_row, "notes"),
    }
    row["proteome_id"] = compute_proteome_id(row)
    return row


def phenotype_labels_for_row(row: Dict[str, object]) -> List[str]:
    labels: List[str] = []
    acid_classes = set(split_multi_value(str(row.get("acidophile_class", ""))))
    alkali_classes = set(split_multi_value(str(row.get("alkaliphile_class", ""))))
    halophile_classes = split_multi_value(str(row.get("halophile_class", "")))

    if "acidotolerant" in acid_classes:
        labels.append("Acidotolerant")
    if str(row.get("acidophile", "")) == "yes":
        if "moderate_acidophile" in acid_classes:
            labels.append("Moderate Acidophile")
        else:
            labels.append("Acidophile")
        if "hyperacidophile" in acid_classes:
            labels.append("Hyperacidophile")
    if str(row.get("acidophilic_alga", "")) == "yes":
        labels.append("Acidophilic Alga")

    if str(row.get("halophile", "")) == "yes":
        if halophile_classes:
            labels.extend(f"Halophile ({title_case_label(value)})" for value in halophile_classes)
        else:
            labels.append("Halophile")
    if str(row.get("saline_adapted_alga", "")) == "yes":
        labels.append("Saline-Adapted Alga")
    if str(row.get("halophilic_alga", "")) == "yes":
        labels.append("Halophilic Alga")

    if str(row.get("alkaliphile", "")) == "yes":
        if "extreme_alkaliphile" in alkali_classes:
            labels.append("Extreme Alkaliphile")
        else:
            labels.append("Alkaliphile")
    if "alkalitolerant" in alkali_classes:
        labels.append("Alkalitolerant")

    if str(row.get("psychrophile", "")) == "yes":
        labels.append("Psychrophile")
    if str(row.get("psychrophilic_alga", "")) == "yes":
        labels.append("Psychrophilic Alga")
    if str(row.get("hyperthermophile", "")) == "yes":
        labels.append("Hyperthermophile")
    elif str(row.get("thermophile", "")) == "yes":
        labels.append("Thermophile")
    if str(row.get("high_temperature_alga", "")) == "yes":
        labels.append("High-Temperature Alga")
    if str(row.get("thermophilic_alga", "")) == "yes":
        labels.append("Thermophilic Alga")
    if str(row.get("thermoacidophilic_alga", "")) == "yes":
        labels.append("Thermoacidophilic Alga")

    ordered = []
    for label in labels:
        if label not in ordered:
            ordered.append(label)
    return ordered


# 5-class labels: subtypes collapse into their broad class (acidotolerant and moderate acidophile -> acidophile,
# alkalitolerant -> alkaliphile); more than one focus = multi-label organism.
def phenotype_focuses_for_row(row: Dict[str, object]) -> List[str]:
    focuses: List[str] = []
    acid_classes = set(split_multi_value(str(row.get("acidophile_class", ""))))
    alkali_classes = set(split_multi_value(str(row.get("alkaliphile_class", ""))))

    if str(row.get("psychrophile", "")) == "yes":
        focuses.append("psychrophile")
    if str(row.get("psychrophilic_alga", "")) == "yes":
        focuses.append("psychrophile")
    if str(row.get("thermophile", "")) == "yes":
        focuses.append("thermophile")
    if str(row.get("halophile", "")) == "yes":
        focuses.append("halophile")
    if (
        str(row.get("high_temperature_alga", "")) == "yes"
        or str(row.get("thermophilic_alga", "")) == "yes"
        or str(row.get("thermoacidophilic_alga", "")) == "yes"
    ):
        focuses.append("thermophile")
    if (
        str(row.get("saline_adapted_alga", "")) == "yes"
        or str(row.get("halophilic_alga", "")) == "yes"
    ):
        focuses.append("halophile")
    if (
        "acidotolerant" in acid_classes
        or "moderate_acidophile" in acid_classes
        or str(row.get("acidophile", "")) == "yes"
        or str(row.get("acidophilic_alga", "")) == "yes"
        or str(row.get("thermoacidophilic_alga", "")) == "yes"
    ):
        focuses.append("acidophile")
    if "alkalitolerant" in alkali_classes or str(row.get("alkaliphile", "")) == "yes":
        focuses.append("alkaliphile")
    return focuses


def aggregate_membership_group(proteome_id: str, rows: List[Dict[str, object]]) -> Dict[str, object]:
    sorted_rows = sort_membership_rows(rows)
    primary = sorted_rows[0]

    halophile_class_values = unique_non_empty(rows, "halophile_class")
    acidophile_class_values = unique_non_empty(rows, "acidophile_class")
    alkaliphile_class_values = unique_non_empty(rows, "alkaliphile_class")
    local_fasta_file, local_fasta_path = choose_canonical_local_path(rows)
    phenotype_focus = join_multi_values(
        focus for row in sorted_rows for focus in phenotype_focuses_for_row(row)
    )
    phenotype_labels = join_multi_values(
        label for row in sorted_rows for label in phenotype_labels_for_row(row)
    )

    aggregated = {
        "proteome_id": proteome_id,
        "dataset_branch": join_multi_values(unique_non_empty(rows, "dataset_branch")),
        "dataset_layer": join_multi_values(unique_non_empty(rows, "dataset_layer")),
        "domain_scope": join_multi_values(unique_non_empty(rows, "domain_scope")),
        "phenotype_focus": phenotype_focus,
        "phenotype_labels": phenotype_labels,
        "branch_membership_count": len(rows),
        "taxon_name": primary.get("taxon_name", ""),
        "ncbi_taxid": first_non_empty(sorted_rows, "ncbi_taxid"),
        "strain": primary.get("strain", ""),
        "strain_aliases": join_multi_values(unique_non_empty(rows, "strain_aliases")),
        "identity_level": join_multi_values(unique_non_empty(rows, "identity_level")),
        "collection_source": join_multi_values(unique_non_empty(rows, "collection_source")),
        "collection_match": join_multi_values(unique_non_empty(rows, "collection_match")),
        "psychrophile": "yes" if any(str(row.get("psychrophile", "")) == "yes" for row in rows) else "no",
        "psychrophilic_alga": "yes"
        if any(str(row.get("psychrophilic_alga", "")) == "yes" for row in rows)
        else "no",
        "thermophile": "yes" if any(str(row.get("thermophile", "")) == "yes" for row in rows) else "no",
        "hyperthermophile": "yes" if any(str(row.get("hyperthermophile", "")) == "yes" for row in rows) else "no",
        "high_temperature_alga": "yes"
        if any(str(row.get("high_temperature_alga", "")) == "yes" for row in rows)
        else "no",
        "thermophilic_alga": "yes"
        if any(str(row.get("thermophilic_alga", "")) == "yes" for row in rows)
        else "no",
        "thermoacidophilic_alga": "yes"
        if any(str(row.get("thermoacidophilic_alga", "")) == "yes" for row in rows)
        else "no",
        "halophile": "yes" if any(str(row.get("halophile", "")) == "yes" for row in rows) else "no",
        "saline_adapted_alga": "yes"
        if any(str(row.get("saline_adapted_alga", "")) == "yes" for row in rows)
        else "no",
        "halophilic_alga": "yes"
        if any(str(row.get("halophilic_alga", "")) == "yes" for row in rows)
        else "no",
        "halophile_class": join_multi_values(halophile_class_values),
        "acidophile": "yes" if any(str(row.get("acidophile", "")) == "yes" for row in rows) else "no",
        "acidophilic_alga": "yes"
        if any(str(row.get("acidophilic_alga", "")) == "yes" for row in rows)
        else "no",
        "acidophile_class": join_multi_values(acidophile_class_values),
        "alkaliphile": "yes" if any(str(row.get("alkaliphile", "")) == "yes" for row in rows) else "no",
        "alkaliphile_class": join_multi_values(alkaliphile_class_values),
        "ogt_celsius": first_non_empty(sorted_rows, "ogt_celsius"),
        "nacl_optimum_M": first_non_empty(sorted_rows, "nacl_optimum_M"),
        "ph_optimum": first_non_empty(sorted_rows, "ph_optimum"),
        "evidence_type": join_multi_values(unique_non_empty(rows, "evidence_type")),
        "evidence_source": join_text_blobs(unique_non_empty(rows, "evidence_source")),
        "include_in_core": "yes" if any(str(row.get("include_in_core", "")) == "yes" for row in rows) else "no",
        "proteome_status": choose_best_status(sorted_rows),
        "proteome_in_manifest": "yes" if any(str(row.get("proteome_in_manifest", "")) == "yes" for row in rows) else "no",
        "protein_count": first_non_empty(sorted_rows, "protein_count"),
        "sequence_source_db": join_multi_values(unique_non_empty(rows, "sequence_source_db")),
        "sequence_source_accession": choose_sequence_accession(sorted_rows),
        "sequence_level": join_multi_values(unique_non_empty(rows, "sequence_level")),
        "retrieval_status": join_multi_values(unique_non_empty(rows, "retrieval_status")),
        "retrieval_notes": join_text_blobs(unique_non_empty(rows, "retrieval_notes")),
        "source_table_file": join_multi_values(unique_non_empty(rows, "source_table_file")),
        "retrieval_tracker_file": join_multi_values(unique_non_empty(rows, "retrieval_tracker_file")),
        "manifest_file": join_multi_values(unique_non_empty(rows, "manifest_file")),
        "merged_fasta_file": join_multi_values(unique_non_empty(rows, "merged_fasta_file")),
        "local_fasta_file": local_fasta_file,
        "local_fasta_path": local_fasta_path,
        "local_fasta_paths": join_multi_values(unique_non_empty(rows, "local_fasta_path")),
        "notes": join_text_blobs(unique_non_empty(rows, "notes")),
    }
    return aggregated


def build_summary(master_rows: List[Dict[str, object]], membership_rows: List[Dict[str, object]]) -> List[Dict[str, object]]:
    grouped: Dict[Tuple[str, str, str], Dict[str, int]] = defaultdict(
        lambda: {"row_count": 0, "protein_count_sum": 0}
    )
    seen_keys: set[Tuple[str, str, str]] = set()

    master_by_id = {str(row["proteome_id"]): row for row in master_rows}
    for membership in membership_rows:
        proteome_id = str(membership["proteome_id"])
        branch = str(membership["dataset_branch"])
        layer = str(membership["dataset_layer"])
        status = str(master_by_id.get(proteome_id, {}).get("proteome_status", membership.get("proteome_status", "")))
        dedupe_key = (proteome_id, branch, status)
        if dedupe_key in seen_keys:
            continue
        seen_keys.add(dedupe_key)

        key = (branch, layer, status)
        grouped[key]["row_count"] += 1
        protein_count = str(master_by_id.get(proteome_id, {}).get("protein_count", "")).strip()
        if protein_count:
            grouped[key]["protein_count_sum"] += int(protein_count)

    summary_rows = []
    for key in sorted(grouped):
        dataset_branch, dataset_layer, proteome_status = key
        values = grouped[key]
        summary_rows.append(
            {
                "dataset_branch": dataset_branch,
                "dataset_layer": dataset_layer,
                "proteome_status": proteome_status,
                "row_count": values["row_count"],
                "protein_count_sum": values["protein_count_sum"],
            }
        )
    return summary_rows


def build_membership_index() -> List[Dict[str, object]]:
    rows: List[Dict[str, object]] = []

    for dataset in DATASETS:
        manifest_path, merged_path = resolve_manifest_and_merged_paths(dataset)
        dataset_with_paths = dict(dataset, manifest=manifest_path, merged_fasta=merged_path)
        curation_rows = read_tsv(dataset["source_table"])
        tracker_rows = read_tsv(dataset["retrieval_tracker"]) if dataset["retrieval_tracker"] else []
        manifest_rows = read_tsv(manifest_path) if manifest_path and manifest_path.exists() else []

        tracker_by_key = build_lookup(tracker_rows)
        tracker_by_relaxed_key = build_relaxed_lookup(tracker_rows)
        manifest_by_key = build_lookup(manifest_rows)
        manifest_by_relaxed_key = build_relaxed_lookup(manifest_rows)
        manifest_by_accession = {
            strip_accession_prefix(row["assembly_accession"]): row for row in manifest_rows
        }
        curation_by_key = build_lookup(curation_rows)
        curation_by_relaxed_key = build_relaxed_lookup(curation_rows)
        curation_by_taxon = build_unique_taxon_lookup(curation_rows)
        tracker_by_taxon = build_unique_taxon_lookup(tracker_rows)

        source_pairs: List[Tuple[Optional[Dict[str, str]], Optional[Dict[str, str]]]] = [
            (row, None) for row in curation_rows
        ]
        seen_source_keys = {relaxed_taxon_key(row) for row in curation_rows}
        for tracker_row in tracker_rows:
            key = relaxed_taxon_key(tracker_row)
            if key in seen_source_keys:
                continue
            curation_taxon_row = curation_by_taxon.get(normalized_taxon_name(tracker_row))
            if curation_taxon_row and not curation_taxon_row.get("strain", "").strip():
                continue
            source_pairs.append((None, tracker_row))
            seen_source_keys.add(key)

        for curation_row, tracker_seed_row in source_pairs:
            source_for_key = curation_row or tracker_seed_row
            if source_for_key is None:
                continue
            exact_key = build_taxon_key(source_for_key)
            relaxed_key = relaxed_taxon_key(source_for_key)

            tracker_row = tracker_seed_row or tracker_by_key.get(exact_key)
            if tracker_row is None:
                tracker_row = tracker_by_relaxed_key.get(relaxed_key)
            if tracker_row is None and curation_row is not None and not curation_row.get("strain", "").strip():
                tracker_row = tracker_by_taxon.get(normalized_taxon_name(curation_row))
            if curation_row is None and tracker_row is not None:
                curation_row = curation_by_key.get(build_taxon_key(tracker_row))
                if curation_row is None:
                    curation_row = curation_by_relaxed_key.get(relaxed_key)
                if curation_row is None and not tracker_row.get("strain", "").strip():
                    curation_row = curation_by_taxon.get(normalized_taxon_name(tracker_row))

            manifest_row = None

            if tracker_row and tracker_row.get("sequence_source_accession"):
                accession = strip_accession_prefix(tracker_row["sequence_source_accession"])
                manifest_row = manifest_by_accession.get(accession)

            if manifest_row is None:
                manifest_row = manifest_by_key.get(exact_key)
            if manifest_row is None:
                manifest_row = manifest_by_relaxed_key.get(relaxed_key)

            rows.append(build_membership_row(dataset_with_paths, curation_row, tracker_row, manifest_row))

    rows.sort(
        key=lambda row: (
            str(row["dataset_layer"]),
            str(row["dataset_branch"]),
            str(row["taxon_name"]),
            str(row["strain"]),
        )
    )
    return rows


def build_master_index(membership_rows: List[Dict[str, object]]) -> List[Dict[str, object]]:
    grouped: Dict[str, List[Dict[str, object]]] = defaultdict(list)
    for row in membership_rows:
        grouped[str(row["proteome_id"])].append(row)

    master_rows = [
        aggregate_membership_group(proteome_id, rows)
        for proteome_id, rows in grouped.items()
    ]
    master_rows.sort(
        key=lambda row: (
            str(row["proteome_status"]) != "downloaded",
            str(row["taxon_name"]),
            str(row["strain"]),
        )
    )
    return master_rows


def print_report(master_rows: List[Dict[str, object]], membership_rows: List[Dict[str, object]]) -> None:
    branch_counter = Counter()
    seen_by_branch: set[Tuple[str, str]] = set()
    for row in membership_rows:
        key = (str(row["proteome_id"]), str(row["dataset_branch"]))
        if key in seen_by_branch:
            continue
        seen_by_branch.add(key)
        branch_counter[str(row["dataset_branch"])] += 1

    status_counter = Counter(str(row["proteome_status"]) for row in master_rows)

    print(f"master_index:     {MASTER_INDEX_PATH}")
    print(f"membership_index:{MEMBERSHIP_INDEX_PATH}")
    print(f"summary:          {SUMMARY_PATH}")
    print(f"rows:             {len(master_rows)}")
    print("rows_by_branch:")
    for branch, count in sorted(branch_counter.items()):
        print(f"  {branch}: {count}")
    print("proteome_status:")
    for status, count in sorted(status_counter.items()):
        print(f"  {status}: {count}")


def main() -> None:
    OUTPUT_DIR.mkdir(exist_ok=True)
    membership_rows = build_membership_index()
    master_rows = build_master_index(membership_rows)
    summary_rows = build_summary(master_rows, membership_rows)

    write_tsv(MEMBERSHIP_INDEX_PATH, membership_rows, MEMBERSHIP_FIELDNAMES)
    write_tsv(MASTER_INDEX_PATH, master_rows, MASTER_FIELDNAMES)
    write_tsv(SUMMARY_PATH, summary_rows, SUMMARY_FIELDNAMES)
    print_report(master_rows, membership_rows)


if __name__ == "__main__":
    main()
