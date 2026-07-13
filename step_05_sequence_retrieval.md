# Step 5: Sequence Retrieval For Accepted Core Rows

Use sequence retrieval only for rows with `include_in_core = yes`.

For the current thermophile pass, that means the `33` accepted rows in the
yes-only export.

## Goal

For each accepted taxon row, attach a reproducible sequence source that can be
used for downstream model work.

Prefer retrieving **proteomes or proteomes linked to the curated strain** rather
than collecting proteins one by one.

## Recommended retrieval unit

Use this priority order:

1. exact strain proteome
2. exact strain genome assembly with translated proteins
3. species proteome only if exact strain retrieval is not available and this is
   clearly documented

## What to add next

For each accepted row, create retrieval metadata such as:

- `sequence_source_db`
- `sequence_source_accession`
- `sequence_level`
- `retrieval_status`
- `retrieval_notes`

Recommended values:

- `sequence_source_db`: `NCBI`, `UniProt`, `GTDB`, `IMG`, `other`
- `sequence_level`: `strain_proteome`, `strain_genome_proteins`,
  `species_proteome`
- `retrieval_status`: `ready`, `review`, `missing`

## Suggested workflow

### 1. Start from the yes-only table

Use the accepted rows only.

Do not retrieve sequences for:

- `review`
- `no`

### 2. Search by NCBI TaxID first

For each accepted row:

- start with `ncbi_taxid`
- prefer an assembly or proteome matching the strain
- if multiple assemblies exist, prefer:
  - reference or representative assemblies
  - complete genomes
  - entries whose strain name matches your curated row

### 3. Record the retrieval decision

For each row, record:

- whether an exact strain proteome was found
- the accession used
- whether the retrieval is exact strain or fallback species level
- anything ambiguous that should stay `review`

### 4. Export FASTA only after the metadata table is stable

First stabilize the retrieval table.

Only then:

- download the selected protein FASTA files
- keep one folder for raw downloads
- keep one merged FASTA for the accepted thermophile set

## Practical first milestone

The next good milestone is:

- a yes-only core table
- a retrieval table with one row per accepted taxon
- exact source accessions for as many rows as possible

After that, sequence download becomes a mechanical step.
