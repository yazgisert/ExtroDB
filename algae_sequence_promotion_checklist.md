# Algae Sequence Promotion Checklist

This note defines when a eukaryotic algae candidate can be promoted from the
paper-screening stage into:

- [core_taxon_table_eukaryotic_algae.tsv](/Users/yazgisert/Documents/Playground/thesis_data_collection/core_taxon_table_eukaryotic_algae.tsv)
- [sequence_retrieval_tracker_eukaryotic_algae_yes.tsv](/Users/yazgisert/Documents/Playground/thesis_data_collection/sequence_retrieval_tracker_eukaryotic_algae_yes.tsv)

The goal is to make the workflow easy to defend in the thesis.

## 1. Short version

Do not move an algae taxon straight from `paper_screening_log` to
`sequence_retrieval_tracker_eukaryotic_algae_yes.tsv`.

Use this order:

1. `AlgaeBase` for accepted-name cleanup
2. `CCAP` or `NCMA` for first-pass identity support
3. primary literature for explicit phenotype evidence
4. `NCBI Taxonomy` for taxon normalization and taxid
5. `NCBI Assembly` and `BioSample` for sequence-linked strain or isolate
   confirmation
6. verify that a downloadable proteome or genome-protein file is actually
   available
7. only then promote to the strict sequence-retrieval tracker

## 2. Why this extra step is needed

`AlgaeBase` is useful for accepted names and synonyms.

It is not enough by itself for sequence retrieval.

To retrieve a sequence responsibly, we need a defensible bridge between:

- the phenotype evidence
- the biological identity
- the actual accession we plan to download

That bridge usually comes from `NCBI Taxonomy`, `NCBI Assembly`, `BioSample`,
and sometimes the primary strain paper.

## 3. Promotion rule

A row can be promoted into
[sequence_retrieval_tracker_eukaryotic_algae_yes.tsv](/Users/yazgisert/Documents/Playground/thesis_data_collection/sequence_retrieval_tracker_eukaryotic_algae_yes.tsv)
only if all of the following are true:

1. the accepted name is stable
2. the phenotype evidence is explicit and citable
3. the taxon passes the algae branch threshold
4. the identity level is stated honestly
5. there is a defensible link to a specific sequence accession
6. the accession leads to a proteome we can actually download and place locally

If one of these is still missing, the row stays in screening or in the core
metadata table with a review note or moves into the algae hold tracker.

## 4. Exact source order

### Step 1: Name cleanup

Use `AlgaeBase` to answer:

- what is the accepted name?
- are there important synonyms?
- should the paper name be normalized before searching other sources?

Output:

- `accepted_name_candidate`
- `accepted_name_source`

### Step 2: First identity pass

Use `CCAP` and `NCMA` to answer:

- is there a strain, isolate, or culture collection record?
- are there aliases that might appear in papers or sequence databases?

Output:

- `strain_candidate`
- `identity_level_candidate`
- `collection_source_checked`
- `collection_match_status`

### Step 3: Phenotype confirmation

Use primary literature to answer:

- is the salinity or temperature value explicit?
- is the value an optimum, best growth condition, or only tolerance?
- does it meet the algae branch threshold?

Output:

- `extracted_nacl_optimum_M`
- `extracted_ogt_celsius`
- `extracted_label`
- `evidence_level_candidate`

### Step 4: Taxonomy confirmation for sequence work

Use `NCBI Taxonomy` to answer:

- what taxid should we use?
- does the accepted name map cleanly to the sequence ecosystem?
- are there alternate names we should watch for?

Output:

- `ncbi_taxid`
- normalized taxon name for retrieval

### Step 5: Sequence-linked identity confirmation

Use `NCBI Assembly` and `BioSample` to answer:

- is there an assembly or proteome for this taxon?
- does the assembly metadata show the same strain or isolate?
- does `BioSample` list strain, isolate, or culture-collection aliases?
- is the record exact-strain, type-material, named study strain, or only
  species-level?

Output:

- `sequence_source_db`
- `sequence_source_accession`
- `sequence_level`
- `retrieval_status`
- `retrieval_notes`
- `download_verified`

## 5. Acceptable sequence-link evidence

Use this order of strength.

### Best case

The phenotype paper and the sequence record point to the same strain or
isolate.

Examples of support:

- same strain code
- same culture-collection identifier
- same isolate designation
- type-material link in assembly metadata

### Good fallback

The phenotype paper uses a named study strain, and `NCBI Assembly` or
`BioSample` lists a matching alias.

This is acceptable if the alias chain is written clearly in the notes.

### Weak fallback

Only a species-level sequence exists.

This can be used only in the exploratory algae branch if:

- the phenotype evidence is strong
- no stable strain mapping is recoverable
- the row is marked honestly as species-level
- the fallback is stated clearly in notes

However, species-level linkage by itself is still not enough for the strict
downloadable tracker unless the actual proteome can be retrieved.

## 6. Promotion decisions

Use these decisions consistently.

### `ready_for_core_table`

Use when:

- accepted name is stable
- phenotype is explicit
- label threshold is met
- identity level is honest

This does not yet require a sequence accession.

### `ready_for_sequence_tracker`

Use when:

- core-table conditions are already satisfied
- `ncbi_taxid` is confirmed
- a specific accession is identified
- sequence linkage is defensible
- a downloadable proteome or genome-protein file is verified

This is the strict `yes` tracker state.

### `hold_below_sequence_tracker`

Use when:

- the core-table conditions are satisfied
- sequence-side evidence exists
- but the download route is not yet verified as a ready usable proteome

### `hold_for_recheck`

Use when:

- numeric phenotype still needs full-text confirmation
- strain mapping is incomplete
- accession exists but identity linkage is weak

### `exclude`

Use when:

- phenotype threshold is not met
- evidence is only habitat inference
- identity is too unstable to explain

## 7. Thesis-safe wording

Simple explanation:

"I first used AlgaeBase to standardize the name, then checked culture
collections for identity support, then used primary literature for the
phenotype, and finally used NCBI Taxonomy plus Assembly and BioSample to
confirm whether a specific sequence accession could be linked defensibly to
that phenotype record."

Shorter version:

"A taxon was promoted to the strict sequence-retrieval tracker only when
phenotype evidence, sequence-linked identity, and a real downloadable proteome
were all defensible."

## 8. Current status of the first three saline-adapted candidates

### Dunaliella tertiolecta

Current status:

- accepted in the exploratory algae core
- explicit salinity evidence is present
- sequence-linked assembly support exists
- current download route is not yet verified as a ready proteome

Decision:

- `hold_below_sequence_tracker`

Needed before promotion:

- verify a real downloadable proteome from the selected accession

### Picochlorum celeri

Current status:

- explicit saline-growth evidence is strong
- first-pass `CCAP` and `NCMA` search did not yield an immediate match

Decision:

- `ready_for_core_table`
- `hold_below_sequence_tracker`

Needed before promotion:

- verify a real downloadable proteome from the current accession

### Picochlorum renovo

Current status:

- accepted in the exploratory algae core
- strain-adjacent sequence support exists through `NREL1`
- the alias caveat remains explicit
- current download route is not yet verified as a ready proteome

Decision:

- `hold_below_sequence_tracker`

Needed before promotion:

- verify a real downloadable proteome from the current accession

## 9. Practical next step

Apply this checklist one taxon at a time in this order:

1. `Picochlorum celeri`
2. `Picochlorum renovo`
3. `Dunaliella tertiolecta`

That order is most efficient because `Picochlorum celeri` currently has the
strongest phenotype support.
