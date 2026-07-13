# Algae Branch Screening Workflow And Templates

## 1. Purpose

This note defines how the eukaryotic algae branch should be restarted in a way
that is:

- compatible with the rest of the thesis workflow
- explicit about what is source-backed versus project-operational
- conservative enough for thesis-grade curation

The algae branch is a separate exploratory branch.
It is not part of the active prokaryotic core.

## 2. The most important starting rule

Do not begin with sequence retrieval.

The correct order is:

1. taxonomic cleanup
2. identity reconciliation
3. phenotype evidence screening
4. accepted metadata table
5. only then sequence retrieval

## 3. Source stack

Use sources in this order.

### 3.1 Taxonomy and synonym cleanup

Use `AlgaeBase` as the first discovery and name-normalization layer.

Use it for:

- accepted names
- synonym cleanup
- taxonomic direction

Do not use it by itself as phenotype proof.

This is a project workflow rule documented in:

- [criteria_sources_for_thesis.md](/Users/yazgisert/Documents/Playground/thesis_data_collection/criteria_sources_for_thesis.md)
- [operational_definitions_extremophile_dataset_curation.md](/Users/yazgisert/Documents/Playground/thesis_data_collection/operational_definitions_extremophile_dataset_curation.md)

### 3.2 Identity anchoring

Use culture collections as the identity layer.

Verified collection sites:

- [CCAP](https://www.ccap.ac.uk/)
- [NCMA](https://ncma.bigelow.org/)

These are useful for:

- strain names
- isolate aliases
- catalogue-level identity support
- later sequence-linking support

### 3.3 Phenotype evidence

Use primary literature first.

Verified example source pages already checked for this branch:

- [Asteromonas gracilis](https://pubmed.ncbi.nlm.nih.gov/16661724/)
- [Dunaliella salina SA 134](https://pubmed.ncbi.nlm.nih.gov/28808229/)
- [Dunaliella tertiolecta](https://pubmed.ncbi.nlm.nih.gov/14610888/)
- [Galdieria sulphuraria](https://pubmed.ncbi.nlm.nih.gov/34202768/)
- [Cyanidioschyzon merolae, context paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC10154280/)
- [Cyanidioschyzon merolae, temperature acclimation paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC5411153/)

Use reviews only as support, not as sole acceptance evidence when a primary
paper should exist.

## 4. What is actually source-backed

The following claims are directly supported by verified source pages.

### 4.1 Asteromonas gracilis

The PubMed abstract states that:

- growth occurs from `0.5 M NaCl` to saturation
- specific growth rate is maximal between `0.5 and 2.5 M NaCl`

Source:

- [PubMed 16661724](https://pubmed.ncbi.nlm.nih.gov/16661724/)

This is a strong example of the type of explicit salinity evidence we want.

### 4.2 Galdieria sulphuraria

The verified PubMed page explicitly describes `Galdieria sulphuraria` as an
`extremophilic` red microalga and states that the paper reports `optimal
conditions`.

Source:

- [PubMed 34202768](https://pubmed.ncbi.nlm.nih.gov/34202768/)

### 4.3 Cyanidioschyzon merolae

The verified PMC context paper explicitly describes `Cyanidioschyzon merolae`
as an `extremophilic red microalga` that grows in `low-pH, high-temperature`
environments.

Source:

- [PMC10154280](https://pmc.ncbi.nlm.nih.gov/articles/PMC10154280/)

The verified PMC temperature paper explicitly treats `C. merolae` as a
thermophilic alga in a temperature-acclimation context.

Source:

- [PMC5411153](https://pmc.ncbi.nlm.nih.gov/articles/PMC5411153/)

### 4.4 What should still be treated cautiously

Some archived algae rows contain project interpretations that should be
re-checked in full text before being accepted again.

Example:

- the current archived `Dunaliella tertiolecta` row should not be treated as
  final numeric optimum evidence until the full paper is re-read carefully

Source page checked:

- [PubMed 14610888](https://pubmed.ncbi.nlm.nih.gov/14610888/)

## 5. What is project-operational rather than universal

The following are explicit project rules.
They are reasonable, but they should be described as operational definitions
in the thesis.

### 5.1 Identity hierarchy

Use this identity preference order:

1. culture-collection strain or isolate
2. named study strain with traceable publication history
3. species-level entry only if no stable strain mapping exists

### 5.2 Habitat words are not enough

Do not accept rows based only on words such as:

- marine
- brackish
- saline environment
- hot spring
- acid spring

### 5.3 Algae salinity labels

Use:

- `saline-adapted alga` for explicit optimum `0.2 <= NaCl equivalent < 0.85 M`
- `halophilic alga` for explicit optimum `>= 0.85 M`

This is a project choice.
It is intentionally stricter for algae than simply calling all marine algae
halophilic.

### 5.4 Algae temperature labels

Use:

- `high-temperature alga` for explicit optimum growth `>= 40 C`
- `thermophilic alga` for explicit optimum growth `>= 45 C`
- `thermoacidophilic alga` for high-temperature growth plus explicit acidic
  growth context

These are operational thesis labels.

## 6. Inclusion criteria for the screening stage

A candidate should move toward acceptance only if all of the following are
true:

- the accepted taxon name is stable enough to defend
- the identity level is clear enough to state honestly
- the phenotype value is explicit in a citable source
- the phenotype passes the algae branch threshold
- the row is not based only on habitat inference

## 7. Decision labels

Use these decisions consistently.

### `include_in_core = yes`

Use `yes` within the algae branch when:

- identity is defensible
- phenotype evidence is explicit
- threshold is satisfied
- source is directly citable

### `include_in_core = review`

Use `review` when:

- the phenotype looks promising
- but the identity mapping is still incomplete
- or the numeric evidence still needs full-text checking

### `include_in_core = no`

Use `no` when:

- the threshold is not met

## 8. Promotion to sequence retrieval

Do not promote directly from the paper-screening log into
[sequence_retrieval_tracker_eukaryotic_algae_yes.tsv](/Users/yazgisert/Documents/Playground/thesis_data_collection/sequence_retrieval_tracker_eukaryotic_algae_yes.tsv).

Use this staged order instead:

1. paper screening
2. core metadata acceptance
3. `NCBI Taxonomy` confirmation
4. `NCBI Assembly` and `BioSample` sequence-link check
5. verify that a downloadable proteome is actually available
6. only then promote into the strict `yes` tracker

The dedicated checklist for this step is documented in:

- [algae_sequence_promotion_checklist.md](/Users/yazgisert/Documents/Playground/thesis_data_collection/algae_sequence_promotion_checklist.md)

If sequence-linked evidence exists but the proteome is not yet actually
downloadable, keep the taxon in:

- [sequence_retrieval_tracker_eukaryotic_algae_hold.tsv](/Users/yazgisert/Documents/Playground/thesis_data_collection/sequence_retrieval_tracker_eukaryotic_algae_hold.tsv)

## 8. How to use the screening logs

The screening logs are first-pass evidence tables.

Each row should capture:

- where the candidate came from
- what name was searched
- what accepted name you currently think is correct
- whether a collection-backed identity exists
- what explicit phenotype value was extracted
- what the next action should be

Use the halophilic and thermophilic screening logs separately so that the
reason for each search remains obvious.

## 9. Recommended first-pass workflow

### Step 1

Pick one focus category:

- `halophilic_algae`
- `thermophilic_algae`

### Step 2

Build a short genus-level candidate list.

Recommended starting genera:

- `Dunaliella`
- `Asteromonas`
- `Galdieria`
- `Cyanidioschyzon`
- `Cyanidium`

### Step 3

For each candidate, do taxonomic cleanup first.

Record:

- searched name
- accepted name
- where the accepted name came from

### Step 4

Check culture collections.

Record:

- whether a collection entry exists
- the collection source
- whether there is a usable strain or isolate match

### Step 5

Read the primary source and extract only explicit values.

Record:

- optimum growth temperature
- optimum salinity
- optimum pH if given
- the phenotype label justified by the source

### Step 6

Assign a follow-up action:

- `collect_identity_support`
- `read_full_text`
- `promote_to_core_table`
- `exclude`

### Step 7

Once a row is stable enough, move it from the screening log into the curated
algae core table.

### Step 8

Only after a stable accepted set exists, build the sequence retrieval layer.

Use two states:

- strict-downloadable rows in [sequence_retrieval_tracker_eukaryotic_algae_yes.tsv](/Users/yazgisert/Documents/Playground/thesis_data_collection/sequence_retrieval_tracker_eukaryotic_algae_yes.tsv)
- sequence-linked but not yet usable rows in [sequence_retrieval_tracker_eukaryotic_algae_hold.tsv](/Users/yazgisert/Documents/Playground/thesis_data_collection/sequence_retrieval_tracker_eukaryotic_algae_hold.tsv)

## 10. Practical first-pass scope

Do not try to build a large algae branch immediately.

A good first milestone is:

- `8-15` screened candidates
- `5-10` strong `review` or `yes` rows
- `2-5` fully identity-stable accepted rows

## 11. Files to use

Screening logs:

- [paper_screening_log_eukaryotic_algae_halophilic.tsv](/Users/yazgisert/Documents/Playground/thesis_data_collection/paper_screening_log_eukaryotic_algae_halophilic.tsv)
- [paper_screening_log_eukaryotic_algae_thermophilic.tsv](/Users/yazgisert/Documents/Playground/thesis_data_collection/paper_screening_log_eukaryotic_algae_thermophilic.tsv)

Curated metadata table:

- [core_taxon_table_eukaryotic_algae.tsv](/Users/yazgisert/Documents/Playground/thesis_data_collection/core_taxon_table_eukaryotic_algae.tsv)

Retrieval tracker:

- [sequence_retrieval_tracker_eukaryotic_algae_yes.tsv](/Users/yazgisert/Documents/Playground/thesis_data_collection/sequence_retrieval_tracker_eukaryotic_algae_yes.tsv)
- [sequence_retrieval_tracker_eukaryotic_algae_hold.tsv](/Users/yazgisert/Documents/Playground/thesis_data_collection/sequence_retrieval_tracker_eukaryotic_algae_hold.tsv)
