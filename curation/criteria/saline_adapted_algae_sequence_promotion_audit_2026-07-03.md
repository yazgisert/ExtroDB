# Saline-Adapted Algae Sequence Promotion Audit

Date: `2026-07-03`

Cleanup note:

This document records the first promotion pass only.

After a stricter cleanup on `2026-07-03`, the project rule for
[sequence_retrieval_tracker_eukaryotic_algae_yes.tsv](/Users/yazgisert/Documents/Playground/thesis_data_collection/sequence_retrieval_tracker_eukaryotic_algae_yes.tsv)
was tightened so that `yes` now means an actually verified downloadable
proteome route.

Rows with sequence-linked evidence but no verified downloadable proteome were
moved into:

- [sequence_retrieval_tracker_eukaryotic_algae_hold.tsv](/Users/yazgisert/Documents/Playground/thesis_data_collection/sequence_retrieval_tracker_eukaryotic_algae_hold.tsv)

This note records the first `NCBI` sequence-side promotion check for the three
current saline-adapted algae candidates.

It follows the workflow documented in:

- [algae_sequence_promotion_checklist.md](/Users/yazgisert/Documents/Playground/thesis_data_collection/algae_sequence_promotion_checklist.md)

## 1. Why this audit exists

The screening log is useful for phenotype-first discovery.

However, promotion into
[sequence_retrieval_tracker_eukaryotic_algae_yes.tsv](/Users/yazgisert/Documents/Playground/thesis_data_collection/sequence_retrieval_tracker_eukaryotic_algae_yes.tsv)
requires more than a promising phenotype paper.

For each candidate, this audit asks:

1. does `NCBI Taxonomy` recognize the taxon cleanly?
2. is there a matching `NCBI Assembly`?
3. does `BioSample` provide strain or isolate support?
4. is the sequence link strong enough for promotion?

## 2. Candidate outcomes

### 2.1 Picochlorum celeri

Sequence-side findings:

- `NCBI Taxonomy` returns taxid `2695429`
- `NCBI` scientific name is `Picochlorum sp. 'celeri'`
- `NCBI Assembly` returns `GCA_010909725.1`
- linked `BioSample` is `SAMN13721886`
- the `BioSample` sample name is `Picochlorum celeri`
- `BioSample` lists `strain = not collected` and `isolate = not collected`

Interpretation:

- phenotype evidence is still strong
- sequence-side support exists
- but the identity is species-level only
- name normalization should explicitly note that the sequence ecosystem uses
  `Picochlorum sp. 'celeri'`

Decision:

- not yet ready for the `yes` sequence tracker
- potentially promotable to the core algae metadata table after name wording is
  normalized

Why it is not yet in sequence retrieval:

- there is no strain-level bridge
- promotion would currently rely on a species-level genome record

### 2.2 Picochlorum renovo

Sequence-side findings:

- `NCBI Taxonomy` returns taxid `2601844`
- `NCBI` scientific name is `Picochlorum renovo (nom. nud.)`
- `NCBI Assembly` returns `GCA_009650465.1`
- linked `BioSample` is `SAMN12509183`
- `BioSample` lists `strain = NREL1`
- the `BioSample` comment states that the alga was isolated and characterized
  at `NREL`

Interpretation:

- sequence-side support clearly exists
- but the current screening row is built around strain candidate `NREL 39-A8`
- the sequence-side record instead points to `NREL1`

Decision after user-directed exception on `2026-07-03`:

- promoted into the exploratory `yes` sequence tracker with an explicit alias
  caveat

Why this was allowed:

- the project has prior exploratory-branch precedent for carrying forward
  defensible strain-adjacent cases
- `NCBI` still provides a concrete taxid, assembly accession, and BioSample
  record for the sequence side
- the mismatch between phenotype-side `NREL 39-A8` and sequence-side `NREL1`
  is now written directly into the retrieval notes rather than hidden

### 2.3 Dunaliella tertiolecta

Sequence-side findings:

- `NCBI Taxonomy` returns taxid `3047`
- `NCBI Assembly` search returns at least two assemblies:
- `GCA_046862625.1`, with `BioSample` support for cultivar `UTEX 999`
- `GCA_031753455.1`, with `BioSource` support for strain `CCAP 19/42`
- `GCA_031753455.1` is labeled `reference genome` in the assembly summary

Interpretation:

- sequence-side support is better than the first quick collection pass
- this taxon does have sequence-linked identity options
- the PubMed abstract itself explicitly confirms the optimum salinity range of
  `0.2-0.5 mol/L NaCl`

Decision after recheck on `2026-07-03`:

- promoted into the exploratory `yes` sequence tracker as a species-level
  fallback case

Why this was allowed:

- the phenotype evidence is explicit in the PubMed abstract itself
- `NCBI Taxonomy` provides a stable taxid
- `NCBI Assembly` provides concrete sequence-linked identities, including
  `GCA_031753455.1` with strain `CCAP 19/42`
- this matches prior exploratory-branch practice of allowing transparent
  species-level fallback retrieval when exact strain linkage is not recoverable

## 3. Current promotion status summary

| Taxon | Taxonomy status | Assembly found | Strain-level bridge | Promotion status |
| --- | --- | --- | --- | --- |
| `Picochlorum celeri` | yes, but normalized as `Picochlorum sp. 'celeri'` | yes | no | moved to hold tracker after strict-download cleanup |
| `Picochlorum renovo` | yes, but marked `nom. nud.` | yes | accepted as alias-adjacent (`NREL1` vs `NREL 39-A8`) | moved to hold tracker after strict-download cleanup |
| `Dunaliella tertiolecta` | yes | yes | species-level fallback accepted | moved to hold tracker after strict-download cleanup |

## 4. Practical conclusion

This check confirms that the promotion workflow is correct:

- `AlgaeBase` should remain the first name-cleanup step
- `CCAP` and `NCMA` remain good first identity checks
- `NCBI Taxonomy`, `Assembly`, and `BioSample` are the right next layer for
  promotion decisions

After the strict-download cleanup on `2026-07-03`, `Picochlorum renovo` no
longer remains in the `yes` tracker.

`Dunaliella tertiolecta` no longer needs to remain phenotype-provisional,
because the PubMed abstract itself explicitly states the `0.2-0.5 mol/L NaCl`
optimum range.

However, it now remains in the hold tracker until a real downloadable proteome
route is verified.

`Picochlorum celeri` remains the strongest unresolved species-level sequence
fallback candidate, but it also remains below the strict downloadable tracker.
