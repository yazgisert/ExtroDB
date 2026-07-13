# AFESM Candidate Table Schema

This schema is for the **candidate-discovery table**, not the final curated
extremophile dataset.

The goal is to preserve AFESM-derived context first, then enrich and filter it
 into a phenotype-backed thesis dataset later.

## Why this table exists

AFESM is most useful for:
- large-scale candidate discovery
- taxonomy context
- biome context
- prioritizing underexplored metagenomic protein clusters

AFESM is not, by itself, the final source of truth for:
- optimal growth temperature
- optimal pH
- optimal salinity
- curated extremophile class labels

So this table should store both:
- what AFESM tells us directly
- what we later infer or add during phenotype curation

## Recommended columns

### Identification

- `afesm_rep_id`
  AFESM representative cluster ID. Primary join key across AFESM-derived files.

- `afesm_member_id`
  Optional member/sequence ID if you later expand from representative-only to
  all-member analysis.

- `source_dataset`
  One of: `AFESM`, `AFDB`, `ESMatlas`, `unknown`.

- `candidate_source`
  Human-readable provenance such as `Yeo_2026_AFESM`.

### Taxonomy context

- `predicted_tax_id`
  Taxonomy ID assigned by AFESM/MMseqs taxonomy workflow.

- `predicted_tax_name`
  Taxonomic name associated with the predicted tax ID.

- `predicted_tax_rank`
  Rank such as `species`, `genus`, `family`, `order`, etc.

- `taxonomy_group_name`
  Higher-level grouping from the AFESM grouping table, if available.

- `taxonomy_specificity_tier`
  Suggested derived field:
  - `high` for species/genus
  - `medium` for family/order
  - `low` for class/phylum/root/unknown

### Biome context

- `biome_lca`
  Lowest common biome assignment from the AFESM biome workflow.

- `n_biomes`
  Number of biome labels contributing to the LCA.

- `biome_lca_depth`
  The LCA depth / specificity length from AFESM.

- `extreme_biome_flag`
  Suggested derived field. `yes` if the biome string suggests potentially
  relevant extreme context, otherwise `no`.

- `extreme_biome_category`
  Suggested normalized category such as:
  - `salt`
  - `marine`
  - `deep_marine_candidate`
  - `host`
  - `soil`
  - `engineered`
  - `unknown`

### Thesis curation fields

- `organism_name_curated`
  Curated organism name after downstream enrichment. May differ from predicted
  AFESM label or remain empty initially.

- `ncbi_taxid_curated`
  Curated stable NCBI TaxID after manual/database-backed reconciliation.

- `extremophile_candidate_type`
  Working hypothesis for follow-up curation, e.g.:
  - `halophile_candidate`
  - `marine_candidate`
  - `thermophile_candidate`
  - `unclassified`

- `evidence_type`
  What currently supports inclusion. Suggested values:
  - `afesm_taxonomy_biome_only`
  - `database_annotation`
  - `literature_supported`
  - `direct_phenotype`

- `growth_temperature_c`
  Optional curated phenotype value.

- `optimal_pH`
  Optional curated phenotype value.

- `optimal_salinity`
  Optional curated phenotype value.

- `polyextremophile_notes`
  Free-text note if the taxon likely belongs to more than one extremophile class.

- `include_in_core_dataset`
  Final decision field:
  - `yes`
  - `no`
  - `review`

- `exclude_reason`
  Why a candidate was excluded, e.g.:
  - `taxonomy_too_unspecific`
  - `no_phenotype_support`
  - `biome_not_relevant`
  - `ambiguous_mapping`

- `review_notes`
  Free-text comment for curation steps, literature links, or unresolved issues.

## Minimum useful export for a first AFESM candidate table

If you want the smallest workable first version, keep these columns:

- `afesm_rep_id`
- `predicted_tax_id`
- `predicted_tax_name`
- `predicted_tax_rank`
- `taxonomy_group_name`
- `biome_lca`
- `n_biomes`
- `biome_lca_depth`
- `extreme_biome_category`
- `extremophile_candidate_type`
- `evidence_type`
- `include_in_core_dataset`
- `review_notes`

## Suggested first-pass interpretation rules

- If `biome_lca` contains `Salt`, mark `extreme_biome_category = salt`.
- If `biome_lca` contains `Marine` or `Aquatic:Marine`, mark it as
  `marine`, not automatically `extremophile`.
- If `predicted_tax_rank` is broader than `family`, do not include it directly
  in the core dataset without additional support.
- If a candidate has only AFESM-derived evidence and no phenotype support, keep
  it in the candidate table but set `include_in_core_dataset = review`.

## Important distinction

This table is for **candidate discovery**.

Your thesis core dataset should be a separate derived table with stronger fields:
- stable taxon identity
- phenotype support
- explicit inclusion criteria
- downstream labels suitable for validation and modeling
