# Step 8: Prokaryotic Core Expansion

The current searchable dataset remains the finalized prokaryotic core.

That means:

- thermophile core remains frozen at the current curated set
- halophile core remains frozen at the original balanced BacDive-backed set
- new expansion work should be generated in separate files first

## Why the next expansion starts with halophiles

Thermophile expansion from the current local TEMPURA plus BacDive workflow is
nearly exhausted.

Halophile expansion is different:

- the original BacDive halophile script intentionally selected a balanced set
  of `11` slight, `11` moderate, and `11` extreme halophiles
- this means the first halophile core was a curated starter core, not the full
  BacDive-supported pool

So the safest next move is:

- keep the same BacDive evidence standard
- continue selecting additional clear-optimum type strains
- store them as a second-pass batch before merging anything

## Current expansion workflow

### 1. Biological expansion batch

Use:

- `expand_halophile_bacdive.py`

This script:

- reuses the original BacDive NaCl parsing and class logic
- excludes the already selected halophile core BacDive IDs
- writes a new second-pass batch instead of overwriting the finalized core

### 2. Retrieval tracker seed

Use:

- `initialize_sequence_tracker.py`

This creates a retrieval tracker from a core-style TSV so the new batch can go
through the same GTDB-first retrieval workflow.

### 3. GTDB-first sequence linking

Use:

- `update_sequence_tracker_from_gtdb.py`

This gives a first pass on which new accepted rows are:

- exact strain-ready
- plausible but review-level
- missing

## Current round 2 outputs

The current second-pass halophile expansion files are:

- `paper_screening_log_halophile_bacdive_expansion_round2.tsv`
- `core_taxon_table_halophile_bacdive_expansion_round2.tsv`
- `bacdive_halophile_shortlist_expansion_round2.tsv`
- `sequence_retrieval_tracker_halophile_bacdive_expansion_round2.tsv`

## Current round 2 status

Biological acceptance:

- `15` new BacDive-backed type strains selected
- balanced as `5` slight, `5` moderate, and `5` extreme halophiles

GTDB retrieval status after alias-aware matching:

- `15` ready
- `0` review
- `0` missing

## Important interpretation

These round 2 rows are not automatically merged into the finalized halophile
core yet.

They should be treated as:

- accepted biological candidates under the same phenotype logic
- sequence-linked expansion rows that are ready for FASTA download and manifest
  building

## Next recommended action

The round 2 batch can now move toward:

- FASTA download
- manifest creation
- merged FASTA creation
- eventual integration into the searchable prokaryotic resource
