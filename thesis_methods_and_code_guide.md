# Thesis Methods And Code Guide

Last updated: 2026-07-05

## 1. Purpose of this guide

This file is a practical explanation of how the thesis dataset was built,
curated, merged, indexed, and prepared for downstream use.

It is meant for writing the thesis methods section and for helping you remember
how the code and data flow actually work.

## 2. Big-picture workflow

The full workflow has been:

1. define phenotype rules
2. screen candidates
3. validate phenotype evidence
4. curate strain-level core tables
5. build sequence retrieval trackers
6. download or place raw protein FASTA files
7. generate manifests and merged FASTAs
8. build the searchable master resource
9. use the deduplicated master index for ML-facing work

## 3. Core data objects in the project

### 3.1 Screening logs

These are the first-pass evidence collection tables.

Examples:

- [paper_screening_log_psychrophile.tsv](/Users/yazgisert/Documents/Playground/thesis_data_collection/paper_screening_log_psychrophile.tsv)
- [paper_screening_log_tempura_bacdive.tsv](/Users/yazgisert/Documents/Playground/thesis_data_collection/paper_screening_log_tempura_bacdive.tsv)
- [paper_screening_log_halophile_bacdive.tsv](/Users/yazgisert/Documents/Playground/thesis_data_collection/paper_screening_log_halophile_bacdive.tsv)
- [paper_screening_log_acidophile_bacdive.tsv](/Users/yazgisert/Documents/Playground/thesis_data_collection/paper_screening_log_acidophile_bacdive.tsv)
- [paper_screening_log_moderate_acidophile_bacdive.tsv](/Users/yazgisert/Documents/Playground/thesis_data_collection/paper_screening_log_moderate_acidophile_bacdive.tsv)
- [paper_screening_log_acidotolerant_bacdive.tsv](/Users/yazgisert/Documents/Playground/thesis_data_collection/paper_screening_log_acidotolerant_bacdive.tsv)

What they store:

- candidate taxon name
- candidate strain
- evidence source
- extracted phenotype values
- preliminary label
- follow-up decision

### 3.2 Core taxon tables

These are the curated organism-level tables that represent accepted or
documented core rows.

Examples:

- [core_taxon_table_psychrophile.tsv](/Users/yazgisert/Documents/Playground/thesis_data_collection/core_taxon_table_psychrophile.tsv)
- [core_taxon_table_tempura_bacdive_yes.tsv](/Users/yazgisert/Documents/Playground/thesis_data_collection/core_taxon_table_tempura_bacdive_yes.tsv)
- [core_taxon_table_halophile_bacdive.tsv](/Users/yazgisert/Documents/Playground/thesis_data_collection/core_taxon_table_halophile_bacdive.tsv)
- [core_taxon_table_acidophile_bacdive.tsv](/Users/yazgisert/Documents/Playground/thesis_data_collection/core_taxon_table_acidophile_bacdive.tsv)
- [core_taxon_table_moderate_acidophile_bacdive.tsv](/Users/yazgisert/Documents/Playground/thesis_data_collection/core_taxon_table_moderate_acidophile_bacdive.tsv)
- [core_taxon_table_acidotolerant_bacdive.tsv](/Users/yazgisert/Documents/Playground/thesis_data_collection/core_taxon_table_acidotolerant_bacdive.tsv)

Important columns include:

- `taxon_name`
- `ncbi_taxid`
- `strain`
- phenotype flags and subclasses
- optimum values (`ogt_celsius`, `nacl_optimum_M`, `ph_optimum`)
- `evidence_type`
- `evidence_source`
- `include_in_core`
- `notes`

### 3.3 Sequence retrieval trackers

These tables connect curated rows to actual sequence sources.

Examples:

- [sequence_retrieval_tracker_psychrophile_yes.tsv](/Users/yazgisert/Documents/Playground/thesis_data_collection/sequence_retrieval_tracker_psychrophile_yes.tsv)
- [sequence_retrieval_tracker_tempura_bacdive_yes.tsv](/Users/yazgisert/Documents/Playground/thesis_data_collection/sequence_retrieval_tracker_tempura_bacdive_yes.tsv)
- [sequence_retrieval_tracker_halophile_bacdive_yes.tsv](/Users/yazgisert/Documents/Playground/thesis_data_collection/sequence_retrieval_tracker_halophile_bacdive_yes.tsv)
- [sequence_retrieval_tracker_acidophile_bacdive_yes.tsv](/Users/yazgisert/Documents/Playground/thesis_data_collection/sequence_retrieval_tracker_acidophile_bacdive_yes.tsv)
- [sequence_retrieval_tracker_moderate_acidophile_bacdive_yes.tsv](/Users/yazgisert/Documents/Playground/thesis_data_collection/sequence_retrieval_tracker_moderate_acidophile_bacdive_yes.tsv)
- [sequence_retrieval_tracker_acidotolerant_bacdive_yes.tsv](/Users/yazgisert/Documents/Playground/thesis_data_collection/sequence_retrieval_tracker_acidotolerant_bacdive_yes.tsv)
- [sequence_retrieval_tracker_eukaryotic_algae_yes.tsv](/Users/yazgisert/Documents/Playground/thesis_data_collection/sequence_retrieval_tracker_eukaryotic_algae_yes.tsv)
- [sequence_retrieval_tracker_eukaryotic_algae_hold.tsv](/Users/yazgisert/Documents/Playground/thesis_data_collection/sequence_retrieval_tracker_eukaryotic_algae_hold.tsv)

Important columns:

- `sequence_source_db`
- `sequence_source_accession`
- `sequence_level`
- `retrieval_status`
- `retrieval_notes`

Typical retrieval logic:

1. try exact strain mapping
2. prefer GTDB when a strong exact match exists
3. otherwise use NCBI assembly/protein accessions
4. document any fallback clearly

For the exploratory algae branch, use one extra split:

- `yes` means the proteome is actually downloadable for local use
- `hold` means the phenotype and sequence link may be acceptable, but the
  download route is not yet verified

### 3.4 Raw protein FASTA folders

Downloaded protein FASTA files are stored by branch.

Folders:

- [psychrophile_raw_protein_fastas](/Users/yazgisert/Documents/Playground/thesis_data_collection/psychrophile_raw_protein_fastas)
- [thermophile_raw_protein_fastas](/Users/yazgisert/Documents/Playground/thesis_data_collection/thermophile_raw_protein_fastas)
- [halophile_raw_protein_fastas](/Users/yazgisert/Documents/Playground/thesis_data_collection/halophile_raw_protein_fastas)
- [acidophile_raw_protein_fastas](/Users/yazgisert/Documents/Playground/thesis_data_collection/acidophile_raw_protein_fastas)
- [moderate_acidophile_raw_protein_fastas](/Users/yazgisert/Documents/Playground/thesis_data_collection/moderate_acidophile_raw_protein_fastas)
- [acidotolerant_raw_protein_fastas](/Users/yazgisert/Documents/Playground/thesis_data_collection/acidotolerant_raw_protein_fastas)
- [algae_raw_protein_fastas](/Users/yazgisert/Documents/Playground/thesis_data_collection/algae_raw_protein_fastas)

These raw folders are the source for manifest building and merged FASTA
construction.

### 3.5 Manifest tables

Manifest tables are structured summaries of which tracker rows currently have a
local FASTA file and how many protein records that file contains.

Examples:

- [sequences_manifest_psychrophile_28.tsv](/Users/yazgisert/Documents/Playground/thesis_data_collection/psychrophile_manifests/sequences_manifest_psychrophile_28.tsv)
- [sequences_manifest_thermophile_68_of_69.tsv](/Users/yazgisert/Documents/Playground/thesis_data_collection/thermophile_manifests/sequences_manifest_thermophile_68_of_69.tsv)
- [sequences_manifest_halophile_46.tsv](/Users/yazgisert/Documents/Playground/thesis_data_collection/halophile_manifests/sequences_manifest_halophile_46.tsv)
- [sequences_manifest_acidophile_21.tsv](/Users/yazgisert/Documents/Playground/thesis_data_collection/acidophile_manifests/sequences_manifest_acidophile_21.tsv)
- [sequences_manifest_moderate_acidophile_15.tsv](/Users/yazgisert/Documents/Playground/thesis_data_collection/moderate_acidophile_manifests/sequences_manifest_moderate_acidophile_15.tsv)
- [sequences_manifest_acidotolerant_9_of_10.tsv](/Users/yazgisert/Documents/Playground/thesis_data_collection/acidotolerant_manifests/sequences_manifest_acidotolerant_9_of_10.tsv)
- `sequences_manifest_eukaryotic_algae_10.tsv` for the current local algae state

### 3.6 Merged FASTA files

Merged FASTA files are concatenated phenotype-specific proteome sets.

Examples:

- [psychrophile_core_28_merged.faa](/Users/yazgisert/Documents/Playground/thesis_data_collection/psychrophile_merged/psychrophile_core_28_merged.faa)
- [thermophile_core_68_of_69_merged.faa](/Users/yazgisert/Documents/Playground/thesis_data_collection/thermophile_merged/thermophile_core_68_of_69_merged.faa)
- [halophile_core_46_merged.faa](/Users/yazgisert/Documents/Playground/thesis_data_collection/halophile_merged/halophile_core_46_merged.faa)
- [acidophile_core_21_merged.faa](/Users/yazgisert/Documents/Playground/thesis_data_collection/acidophile_merged/acidophile_core_21_merged.faa)
- [moderate_acidophile_core_15_merged.faa](/Users/yazgisert/Documents/Playground/thesis_data_collection/moderate_acidophile_merged/moderate_acidophile_core_15_merged.faa)
- [acidotolerant_core_9_of_10_merged.faa](/Users/yazgisert/Documents/Playground/thesis_data_collection/acidotolerant_merged/acidotolerant_core_9_of_10_merged.faa)
- [acidophile_combined_core_36_merged.faa](/Users/yazgisert/Documents/Playground/thesis_data_collection/acidophile_merged/acidophile_combined_core_36_merged.faa)
- `eukaryotic_algae_core_10_merged.faa` for the current local algae state

## 4. Main phenotype-specific workflows

### 4.1 Thermophiles

Primary workflow:

- candidate discovery from TEMPURA
- phenotype checking and strain validation in BacDive
- sequence tracker population
- GTDB-first retrieval linking
- local FASTA integration

Current scope note:

- this workflow is the hot-side temperature workflow only
- psychrophiles should be handled as a separate psychrophile workflow or under
  a renamed `temperature-related` umbrella
- cold-adapted rows should not be added into tables or folders that still carry
  the explicit `thermophile` name

Relevant files/scripts:

- [curate_tempura_bacdive.py](/Users/yazgisert/Documents/Playground/thesis_data_collection/curate_tempura_bacdive.py)
- [prepare_tempura_second_pass.py](/Users/yazgisert/Documents/Playground/thesis_data_collection/prepare_tempura_second_pass.py)
- [core_taxon_table_tempura_bacdive_yes.tsv](/Users/yazgisert/Documents/Playground/thesis_data_collection/core_taxon_table_tempura_bacdive_yes.tsv)
- [sequence_retrieval_tracker_tempura_bacdive_yes.tsv](/Users/yazgisert/Documents/Playground/thesis_data_collection/sequence_retrieval_tracker_tempura_bacdive_yes.tsv)

### 4.1A Psychrophiles

Current branch state:

- the psychrophile branch now exists as a separate curated branch
- it uses its own screening log, core table, tracker, raw FASTA folder,
  manifest folder, and merged-output folder
- the current psychrophile core now contains `28` accepted rows in the live
  branch files
- the current retrieval tracker marks those `28` rows as `ready` for download

Starter files:

- [paper_screening_log_psychrophile.tsv](/Users/yazgisert/Documents/Playground/thesis_data_collection/paper_screening_log_psychrophile.tsv)
- [core_taxon_table_psychrophile.tsv](/Users/yazgisert/Documents/Playground/thesis_data_collection/core_taxon_table_psychrophile.tsv)
- [sequence_retrieval_tracker_psychrophile_yes.tsv](/Users/yazgisert/Documents/Playground/thesis_data_collection/sequence_retrieval_tracker_psychrophile_yes.tsv)

### 4.2 Halophiles

Primary workflow:

- BacDive-based discovery and validation
- balanced class-aware selection in the earlier core
- later expansion using new BacDive-backed candidates

Relevant files/scripts:

- [curate_halophile_bacdive.py](/Users/yazgisert/Documents/Playground/thesis_data_collection/curate_halophile_bacdive.py)
- [expand_halophile_bacdive.py](/Users/yazgisert/Documents/Playground/thesis_data_collection/expand_halophile_bacdive.py)
- [core_taxon_table_halophile_bacdive.tsv](/Users/yazgisert/Documents/Playground/thesis_data_collection/core_taxon_table_halophile_bacdive.tsv)
- [sequence_retrieval_tracker_halophile_bacdive_yes.tsv](/Users/yazgisert/Documents/Playground/thesis_data_collection/sequence_retrieval_tracker_halophile_bacdive_yes.tsv)

### 4.3 Acid-related branches

Acid-related work was split into three branches:

- strict acidophiles
- moderate acidophiles
- acidotolerants

Hyperacidophiles are stored as a subclass of the strict acidophile branch.

Relevant files:

- [core_taxon_table_acidophile_bacdive.tsv](/Users/yazgisert/Documents/Playground/thesis_data_collection/core_taxon_table_acidophile_bacdive.tsv)
- [core_taxon_table_moderate_acidophile_bacdive.tsv](/Users/yazgisert/Documents/Playground/thesis_data_collection/core_taxon_table_moderate_acidophile_bacdive.tsv)
- [core_taxon_table_acidotolerant_bacdive.tsv](/Users/yazgisert/Documents/Playground/thesis_data_collection/core_taxon_table_acidotolerant_bacdive.tsv)

## 5. Criteria and literature documentation

The criteria and justification layer lives mainly in:

- [criteria_sources_for_thesis.md](/Users/yazgisert/Documents/Playground/thesis_data_collection/criteria_sources_for_thesis.md)
- [operational_definitions_extremophile_dataset_curation.md](/Users/yazgisert/Documents/Playground/thesis_data_collection/operational_definitions_extremophile_dataset_curation.md)
- [extremophile_core_dataset_inclusion_criteria.md](/Users/yazgisert/Documents/Playground/thesis_data_collection/extremophile_core_dataset_inclusion_criteria.md)

When writing the thesis, you can treat these as:

- criteria source notes
- operational definition notes
- minimum-inclusion policy notes

## 6. Sequence and FASTA processing scripts

### 6.1 `build_extremophile_manifests.py`

File:

- [build_extremophile_manifests.py](/Users/yazgisert/Documents/Playground/thesis_data_collection/build_extremophile_manifests.py)

Purpose:

- checks which tracker rows have a matching local FASTA
- counts protein records in each FASTA
- builds manifest tables
- creates merged FASTA files

Key logic:

1. read a sequence retrieval tracker
2. normalize accession names
3. look for matching FASTA files in the correct raw folder
4. count `>` lines as protein records
5. write manifest TSV
6. concatenate FASTAs into a merged output

Typical command:

```bash
python3 /Users/yazgisert/Documents/Playground/thesis_data_collection/build_extremophile_manifests.py
```

### 6.2 `common.py`

File:

- [common.py](/Users/yazgisert/Documents/Playground/thesis_data_collection/common.py)

Purpose:

- shared TSV reading/writing
- accession normalization
- FASTA record counting
- FASTA merging
- path normalization

This file is the helper layer used by multiple scripts.

## 7. Searchable resource and web app

### 7.1 Why it exists

The searchable resource was built so that the curated organism-level dataset
could be:

- queried more easily
- inspected visually
- exported more safely
- connected to local FASTA files and merged outputs

### 7.2 Important change: deduplicated master index

Originally, the searchable layer effectively behaved like a branch-level row
table. That was fine for curation, but not ideal for ML because the same
proteome could appear more than once under different labels.

This was solved by switching to a hybrid model:

- a deduplicated master proteome index
- a separate membership/provenance table

Files:

- [extremophile_master_index.tsv](/Users/yazgisert/Documents/Playground/thesis_data_collection/searchable_resource/extremophile_master_index.tsv)
- [extremophile_membership_index.tsv](/Users/yazgisert/Documents/Playground/thesis_data_collection/searchable_resource/extremophile_membership_index.tsv)
- [extremophile_master_summary.tsv](/Users/yazgisert/Documents/Playground/thesis_data_collection/searchable_resource/extremophile_master_summary.tsv)

### 7.3 `build_searchable_extremophile_index.py`

File:

- [build_searchable_extremophile_index.py](/Users/yazgisert/Documents/Playground/thesis_data_collection/build_searchable_extremophile_index.py)

What it now does:

1. reads all curated phenotype branches
2. links curation rows to retrieval tracker rows
3. links retrieval rows to manifest rows
4. builds a branch-level membership table
5. collapses shared proteomes into one deduplicated proteome row
6. preserves branch and label overlap as multi-value fields
7. writes a summary table

Important idea:

- `proteome_id` is the unique proteome key
- one `proteome_id` can have multiple phenotype labels and multiple source
  branches

This is the reason the master index is safe to use for ML-facing counting.

### 7.4 `searchable_resource_webapp.py`

File:

- [searchable_resource_webapp.py](/Users/yazgisert/Documents/Playground/thesis_data_collection/searchable_resource_webapp.py)

Purpose:

- serves the local app
- exposes `/api/index` and `/api/summary`
- serves local FASTA and manifest files through `/files/...`

Typical command:

```bash
python3 /Users/yazgisert/Documents/Playground/thesis_data_collection/searchable_resource_webapp.py
```

### 7.5 Front-end app

Files:

- [index.html](/Users/yazgisert/Documents/Playground/thesis_data_collection/searchable_resource_app/index.html)
- [app.js](/Users/yazgisert/Documents/Playground/thesis_data_collection/searchable_resource_app/app.js)
- [README.md](/Users/yazgisert/Documents/Playground/thesis_data_collection/searchable_resource_app/README.md)

What the app now supports:

- filtering by branch
- filtering by phenotype focus
- filtering by phenotype subclasses
- showing overlap labels on a single unique proteome row
- linking to local FASTA and merged outputs

Because the app now uses the deduplicated master index, one organism can still
be visible from multiple labels without being counted twice.

## 8. Recommended way to describe the method in the thesis

You can describe the pipeline in these layers:

### Biological curation layer

- literature-backed phenotype definitions
- database-assisted candidate discovery
- strain-aware acceptance rules

### Sequence-linking layer

- accepted rows linked to GTDB or NCBI
- raw protein FASTA files collected locally
- manifest and merged FASTA generation

### Resource-construction layer

- branch-level provenance preserved
- unique proteome rows deduplicated
- searchable local interface built on top

### ML-preparation layer

- master index contains one row per unique proteome
- phenotype overlap stored as annotations, not duplicates

## 9. Current numbers you can cite from the repo

Current deduplicated master resource:

- `145` unique proteomes
- `142` downloaded proteomes
- `446,497` downloaded protein sequences total
- `18` proteomes with overlap across more than one phenotype branch

Current branch-level downloaded sets:

- thermophile: `68 / 69`, `152,071` proteins
- halophile: `46 / 46`, `180,326` proteins
- acidophile: `21 / 21`, `50,761` proteins
- moderate acidophile: `15 / 15`, `55,930` proteins
- acidotolerant: `10 / 10`, `52,214` proteins

## 10. Which files to use for what

Use these files for thesis writing:

- biological criteria: [criteria_sources_for_thesis.md](/Users/yazgisert/Documents/Playground/thesis_data_collection/criteria_sources_for_thesis.md)
- progress summary: [thesis_progress_summary_for_supervisor.md](/Users/yazgisert/Documents/Playground/thesis_data_collection/thesis_progress_summary_for_supervisor.md)
- methods/code memory aid: [thesis_methods_and_code_guide.md](/Users/yazgisert/Documents/Playground/thesis_data_collection/thesis_methods_and_code_guide.md)

Use these files for analysis and export:

- ML-safe master index: [extremophile_master_index.tsv](/Users/yazgisert/Documents/Playground/thesis_data_collection/searchable_resource/extremophile_master_index.tsv)
- branch provenance: [extremophile_membership_index.tsv](/Users/yazgisert/Documents/Playground/thesis_data_collection/searchable_resource/extremophile_membership_index.tsv)
- branch-level sequence assets: manifests and merged FASTAs in the phenotype
  folders

## 11. Important caution

Some older step documents in the repo were written before the acidophile and
acidotolerant branches were integrated into the active resource. They are still
useful historically, but the current ground truth is the present TSVs, merged
FASTAs, and searchable resource outputs.
