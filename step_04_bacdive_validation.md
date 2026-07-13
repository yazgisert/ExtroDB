# Step 4: Validate TEMPURA Candidates In BacDive

Use BacDive as the main **strain-level phenotype validation** source for the
first curation pass.

Do this for the `10` TEMPURA candidates you choose:

- `5` hyperthermophiles
- `5` standard thermophiles

## Goal

For each chosen TEMPURA row, confirm in BacDive:

- the strain identity
- the taxonomic name
- whether a temperature entry exists
- whether the temperature evidence supports the thermophile label
- whether the row should stay `yes`, `review`, or `no`

## BacDive workflow

### 1. Start with the BacDive main search

Search by:

- species name first
- then species name plus strain
- if needed, NCBI TaxID

Useful examples:

- `Thermus thermophilus`
- `Thermus thermophilus HB8`
- `taxid:274`

### 2. Open the most likely strain entry

Prefer the entry that matches:

- the exact species
- the exact strain designation
- the clearest taxonomic record

If there are multiple entries:

- prefer the exact strain match
- prefer a type strain if strain-level certainty is otherwise unclear
- if the strain is still ambiguous, mark the row `review`

### 3. Check these fields in the strain detail page

#### Identity fields

Look under name and taxonomy:

- full scientific name
- strain designation
- species
- genus

#### Temperature fields

Look under culture and growth conditions:

- `culture_temp`
- especially entries where `test_type = optimum`
- also note `growth`, `maximum`, or `minimum` if shown

#### Supporting references

Look under literature/references:

- DOI
- title
- year

## What to record

Update the screening log notes for each checked row with:

- whether exact strain match was found
- whether BacDive has an optimum temperature entry
- whether the BacDive value agrees with TEMPURA
- the BacDive strain page or reference DOI if available

When you move a row into the core table, fill:

- `taxon_name`
- `ncbi_taxid`
- `strain`
- `thermophile`
- `hyperthermophile`
- `ogt_celsius`
- `evidence_type`
- `evidence_source`
- `include_in_core`
- `notes`

## Decision rules

### Mark `yes`

Use `include_in_core = yes` if:

- the BacDive strain identity matches clearly
- BacDive or another strong source gives temperature evidence
- the temperature supports the threshold

Recommended first-pass core evidence types:

- `curated_database`
- `review`
- `experimental`

### Mark `review`

Use `include_in_core = review` if:

- species matches but strain does not match clearly
- BacDive has related temperature evidence but not a clean optimum
- BacDive and TEMPURA disagree
- the value is near the threshold

### Mark `no`

Use `include_in_core = no` if:

- the strain identity appears wrong
- BacDive contradicts the thermophile claim
- no usable support remains after checking

## Practical interpretation for thermophiles

- `45 <= OGT < 80` -> `thermophile = yes`, `hyperthermophile = no`
- `OGT >= 80` -> `thermophile = yes`, `hyperthermophile = yes`

## Keep this simple for the first pass

Do not try to resolve every literature discrepancy immediately.

For the first round:

- confirm identity
- confirm temperature evidence
- make a defendable `yes/review/no` call

That is enough to get your first core rows started.
