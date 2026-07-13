#!/usr/bin/env python3

import csv
import gzip
import io
import re
import sys
import urllib.request
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple


GTDB_RELEASE_BASE = "https://data.gtdb.ecogenomic.org/releases/latest"
METADATA_FILES = [
    f"{GTDB_RELEASE_BASE}/ar53_metadata.tsv.gz",
    f"{GTDB_RELEASE_BASE}/bac120_metadata.tsv.gz",
]


def normalize(text: str) -> str:
    return re.sub(r"[^a-z0-9]", "", text.lower())


def split_variants(value: str) -> List[str]:
    if not value:
        return []
    parts = re.split(r"[;,|]", value)
    variants = [part.strip() for part in parts if part.strip()]
    return variants


def target_strain_variants(target_row: Dict[str, str]) -> List[str]:
    values: List[str] = []
    primary = (target_row.get("strain") or "").strip()
    if primary:
        values.append(primary)
    aliases = target_row.get("strain_aliases", "") or ""
    values.extend(split_variants(aliases))

    seen = set()
    ordered: List[str] = []
    for value in values:
        if value not in seen:
            ordered.append(value)
            seen.add(value)
    return ordered


def field_values(row: Dict[str, str]) -> List[str]:
    values: List[str] = []
    for key in (
        "ncbi_strain_identifiers",
        "ncbi_isolate",
        "ncbi_organism_name",
        "ncbi_assembly_name",
        "ncbi_taxonomy",
        "ncbi_taxonomy_unfiltered",
    ):
        raw = row.get(key, "") or ""
        values.extend(split_variants(raw))
        if raw:
            values.append(raw)
    return [value for value in values if value]


def exact_strain_match(strain: str, values: Iterable[str]) -> bool:
    target = normalize(strain)
    if not target:
        return False
    for value in values:
        if normalize(value) == target:
            return True
    return False


def partial_strain_match(strain: str, values: Iterable[str]) -> bool:
    target = normalize(strain)
    if not target:
        return False
    for value in values:
        norm = normalize(value)
        if not norm:
            continue
        if target in norm or norm in target:
            return True
    return False


def species_name_match(target_species: str, metadata_row: Dict[str, str]) -> bool:
    target = target_species.strip().lower()
    species_tax = metadata_row.get("ncbi_taxonomy", "") or ""
    organism = metadata_row.get("ncbi_organism_name", "") or ""
    return target in species_tax.lower() or organism.lower().startswith(target)


def score_candidate(target_row: Dict[str, str], metadata_row: Dict[str, str]) -> Tuple[int, List[str]]:
    reasons: List[str] = []
    score = 0

    curated_taxid = (target_row.get("ncbi_taxid") or "").strip()
    species_taxid = (metadata_row.get("ncbi_species_taxid") or "").strip()
    ncbi_taxid = (metadata_row.get("ncbi_taxid") or "").strip()
    target_variants = target_strain_variants(target_row)
    values = field_values(metadata_row)

    if curated_taxid and curated_taxid == species_taxid:
        score += 40
        reasons.append("species_taxid")
    elif curated_taxid and curated_taxid == ncbi_taxid:
        score += 30
        reasons.append("ncbi_taxid")

    if species_name_match(target_row["taxon_name"], metadata_row):
        score += 25
        reasons.append("species_name")

    exact_reason = None
    partial_reason = None
    for variant in target_variants:
        if exact_strain_match(variant, values):
            exact_reason = variant
            break
    if exact_reason is not None:
        score += 120
        reasons.append(f"exact_strain:{exact_reason}")
    else:
        for variant in target_variants:
            if partial_strain_match(variant, values):
                partial_reason = variant
                break
        if partial_reason is not None:
            score += 70
            reasons.append(f"partial_strain:{partial_reason}")

    if (metadata_row.get("gtdb_representative") or "").strip().lower() in {"t", "true", "yes"}:
        score += 8
        reasons.append("gtdb_representative")

    type_designation = (metadata_row.get("gtdb_type_designation_ncbi_taxa") or "").lower()
    if "type strain" in type_designation:
        score += 8
        reasons.append("type_strain")

    assembly_level = (metadata_row.get("ncbi_assembly_level") or "").lower()
    if "complete" in assembly_level:
        score += 6
        reasons.append("complete_genome")
    elif "chromosome" in assembly_level:
        score += 4
        reasons.append("chromosome")

    refseq_category = (metadata_row.get("ncbi_refseq_category") or "").lower()
    if "reference genome" in refseq_category:
        score += 5
        reasons.append("reference_genome")
    elif "representative genome" in refseq_category:
        score += 4
        reasons.append("representative_genome")

    return score, reasons


def iter_metadata_rows(url: str) -> Iterable[Dict[str, str]]:
    with urllib.request.urlopen(url) as response:
        with gzip.GzipFile(fileobj=response) as gz:
            text_stream = io.TextIOWrapper(gz, encoding="utf-8", newline="")
            reader = csv.DictReader(text_stream, delimiter="\t")
            for row in reader:
                yield row


def load_targets(path: Path) -> List[Dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def pick_best_candidates(targets: List[Dict[str, str]]) -> Dict[Tuple[str, str], Dict[str, object]]:
    by_key = {(row["taxon_name"], row["strain"]): row for row in targets}
    taxids = {row["ncbi_taxid"] for row in targets if row.get("ncbi_taxid")}
    species_names = {row["taxon_name"].strip().lower() for row in targets}
    candidates: Dict[Tuple[str, str], List[Tuple[int, List[str], Dict[str, str]]]] = defaultdict(list)

    for url in METADATA_FILES:
        for metadata_row in iter_metadata_rows(url):
            row_species_taxid = (metadata_row.get("ncbi_species_taxid") or "").strip()
            row_ncbi_taxid = (metadata_row.get("ncbi_taxid") or "").strip()
            row_species_hit = False

            if row_species_taxid and row_species_taxid in taxids:
                row_species_hit = True
            elif row_ncbi_taxid and row_ncbi_taxid in taxids:
                row_species_hit = True
            else:
                organism = (metadata_row.get("ncbi_organism_name") or "").strip().lower()
                taxonomy = (metadata_row.get("ncbi_taxonomy") or "").strip().lower()
                if any(name in organism or name in taxonomy for name in species_names):
                    row_species_hit = True

            if not row_species_hit:
                continue

            for key, target_row in by_key.items():
                target_name = target_row["taxon_name"].strip().lower()
                if target_name not in (metadata_row.get("ncbi_organism_name") or "").lower() and target_name not in (
                    metadata_row.get("ncbi_taxonomy") or ""
                ).lower():
                    if target_row.get("ncbi_taxid") not in {row_species_taxid, row_ncbi_taxid}:
                        continue

                score, reasons = score_candidate(target_row, metadata_row)
                if score <= 0:
                    continue
                candidates[key].append((score, reasons, metadata_row))

    best: Dict[Tuple[str, str], Dict[str, object]] = {}
    for key, entries in candidates.items():
        entries.sort(
            key=lambda item: (
                item[0],
                (item[2].get("gtdb_representative") or "").strip().lower() in {"t", "true", "yes"},
                item[2].get("ncbi_assembly_level") or "",
            ),
            reverse=True,
        )
        best[key] = {
            "score": entries[0][0],
            "reasons": entries[0][1],
            "row": entries[0][2],
            "candidate_count": len(entries),
        }
    return best


def classify_sequence_level(target_row: Dict[str, str], metadata_row: Dict[str, str]) -> Tuple[str, str]:
    values = field_values(metadata_row)
    for variant in target_strain_variants(target_row):
        if exact_strain_match(variant, values):
            return "strain_genome_proteins", "ready"
    for variant in target_strain_variants(target_row):
        if partial_strain_match(variant, values):
            return "strain_genome_proteins", "review"
    return "species_genome_proteins", "review"


def build_note(target_row: Dict[str, str], metadata_row: Dict[str, str], reasons: List[str], status: str) -> str:
    curated_taxid = (target_row.get("ncbi_taxid") or "").strip()
    gtdb_species_taxid = (metadata_row.get("ncbi_species_taxid") or "").strip()
    gtdb_taxid = (metadata_row.get("ncbi_taxid") or "").strip()
    accession = metadata_row.get("accession") or metadata_row.get("ncbi_genbank_assembly_accession") or ""
    strain_ids = metadata_row.get("ncbi_strain_identifiers") or metadata_row.get("ncbi_isolate") or ""

    note_parts = [f"GTDB accession {accession}"]
    if strain_ids:
        note_parts.append(f"strain identifiers: {strain_ids}")
    if curated_taxid and curated_taxid != gtdb_species_taxid:
        if gtdb_species_taxid:
            note_parts.append(
                f"curated species taxid {curated_taxid} differs from GTDB ncbi_species_taxid {gtdb_species_taxid}"
            )
        elif curated_taxid != gtdb_taxid and gtdb_taxid:
            note_parts.append(f"curated taxid {curated_taxid} differs from GTDB ncbi_taxid {gtdb_taxid}")
    if reasons:
        note_parts.append(f"match basis: {', '.join(reasons)}")
    if status == "review":
        note_parts.append("manual check recommended")
    return "; ".join(note_parts) + "."


def update_tracker(tracker_path: Path, best_matches: Dict[Tuple[str, str], Dict[str, object]]) -> List[Dict[str, str]]:
    with tracker_path.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))

    for row in rows:
        key = (row["taxon_name"], row["strain"])
        match = best_matches.get(key)
        if not match:
            row["sequence_source_db"] = "GTDB"
            row["sequence_source_accession"] = ""
            row["sequence_level"] = ""
            row["retrieval_status"] = "missing"
            row["retrieval_notes"] = "No GTDB genome matched this curated row by species taxid/name."
            continue

        metadata_row = match["row"]  # type: ignore[index]
        reasons = match["reasons"]  # type: ignore[index]
        sequence_level, status = classify_sequence_level(row, metadata_row)  # type: ignore[arg-type]

        row["sequence_source_db"] = "GTDB"
        row["sequence_source_accession"] = metadata_row.get("accession") or metadata_row.get(
            "ncbi_genbank_assembly_accession", ""
        )
        row["sequence_level"] = sequence_level
        row["retrieval_status"] = status
        row["retrieval_notes"] = build_note(row, metadata_row, reasons, status)

    with tracker_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0].keys(), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    return rows


def main() -> int:
    if len(sys.argv) != 3:
        print(
            "Usage: update_sequence_tracker_from_gtdb.py CORE_YES_TSV TRACKER_TSV",
            file=sys.stderr,
        )
        return 1

    core_yes_path = Path(sys.argv[1])
    tracker_path = Path(sys.argv[2])

    targets = load_targets(core_yes_path)
    best_matches = pick_best_candidates(targets)
    updated_rows = update_tracker(tracker_path, best_matches)

    ready = sum(1 for row in updated_rows if row["retrieval_status"] == "ready")
    review = sum(1 for row in updated_rows if row["retrieval_status"] == "review")
    missing = sum(1 for row in updated_rows if row["retrieval_status"] == "missing")
    print(f"ready={ready}\nreview={review}\nmissing={missing}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
