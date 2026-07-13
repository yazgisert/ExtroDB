# Step 1: Build The Core Taxon List

This is the first concrete data-collection step for the thesis.

Do not start by collecting arbitrary protein sequences.

Start with a **curated organism/taxon list** for the first phenotype-backed core
dataset. Sequence collection comes only after this table is stable enough.

## Goal of step 1

Create a first-pass list of taxa that can be labeled with enough confidence as:

- `thermophile`
- `halophile`

These are the safest first categories because they are easier to justify
biologically and usually easier to support with explicit phenotype evidence.

## What to collect first

For each candidate strain or taxon, record:

- taxon name
- strain identifier if available
- NCBI TaxID
- `thermophile` yes/no
- `hyperthermophile` yes/no
- `halophile` yes/no
- `halophile_class` slight/moderate/extreme/null
- `ogt_celsius`
- `nacl_optimum_M`
- `evidence_type`
- `evidence_source`
- `include_in_core`
- notes

Use the TSV template in `core_taxon_table_template.tsv`.

## Recommended evidence standard

### Include now

Use `include_in_core = yes` only if:

- strain or taxon identity is stable
- NCBI TaxID is available or can be reconciled
- the extremophile label is explicit enough to defend
- the evidence source is clear
- evidence is `experimental`, `curated_database`, or `review`
- the value satisfies the operational threshold

### Keep for review

Use `include_in_core = review` if:

- the source suggests the taxon is relevant
- but phenotype support is still incomplete
- evidence is only `genus_inference`
- values are close to the threshold
- strain identity or source agreement is ambiguous

### Exclude for now

Use `include_in_core = no` if:

- taxonomy is too vague
- evidence is too weak
- the mapping between organism and label is ambiguous
- evidence is only `biome_inference`
- the threshold is not satisfied

## Best source strategy

Do not read papers sequence by sequence.

Instead, start from:

- review papers about thermophiles and halophiles
- trusted database annotations
- organism-level phenotype reports

The goal of this step is a **defensible strain-aware taxon table**, not full
coverage.

## Fixed thresholds for the first version

- `thermophile`: `ogt_celsius >= 45`
- `hyperthermophile`: `ogt_celsius >= 80`
- exclude thermotolerant mesophiles with `40 <= OGT < 45`
- `halophile`: `nacl_optimum_M > 0.2`
- `slight halophile`: `0.2 <= NaCl optimum < 0.85`
- `moderate halophile`: `0.85 <= NaCl optimum < 3.4`
- `extreme halophile`: `NaCl optimum >= 3.4`
- exclude halotolerant organisms with salt tolerance but optimum below `0.2 M`

## Evidence vocabulary for the first table

Use these exact `evidence_type` values:

- `experimental`
- `curated_database`
- `review`
- `genus_inference`
- `biome_inference`

## Important unit of curation

Each row should correspond to a **strain** when possible, not just a species.

If only species-level evidence is available:

- prefer the type strain if it can be identified
- keep the row
- note the missing strain-level confirmation in `notes`

## Suggested first milestone

A strong first milestone is:

- 20 to 50 rows total
- mostly strain or species level
- each row documented with evidence
- clear separation between `yes`, `review`, and `no`

Once that exists, move to sequence retrieval for the `yes` set only.
