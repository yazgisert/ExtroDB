# Proposal: Separate Eukaryotes Dataset Branch For The Thesis

## Recommendation

Add a **separate eukaryotic exploratory branch** to the thesis dataset rather
than merging eukaryotes into the active prokaryotic core.

This keeps the current prokaryotic resource methodologically strict and
ML-safe, while still allowing the thesis to demonstrate that the framework can
be extended to a second biological domain with different evidence and identity
constraints.

## Why a separate branch is the right design

The current core dataset is built around:

- explicit phenotype values
- strain-aware identity matching
- traceable proteome retrieval
- branch-aware but deduplicated downstream integration

That logic works well for prokaryotic taxa because culture collections,
BacDive/TEMPURA-style resources, and genome/proteome mappings are often easier
to reconcile at strain level.

For eukaryotes, especially algae, the situation is different:

- strain naming is less consistently recoverable across sources
- species-level physiological studies are more common than strain-ready
  proteome links
- phenotype evidence and sequence identity often come from different source
  layers
- marine tolerance and true halophily are easier to confuse

Because of this, mixing eukaryotes into the current core would weaken the
criteria of the main dataset. A separate branch solves that problem.

## Proposed thesis framing

The thesis should define two dataset layers:

1. **Primary curated prokaryotic core**
   This remains the main thesis dataset and the main ML-facing resource.
2. **Secondary eukaryotic exploratory branch**
   This is a method extension and comparative branch, not the main benchmark.

This framing lets you say:

- the core thesis deliverable remains complete without eukaryotes
- the eukaryote branch demonstrates extensibility of the curation framework
- negative or partial eukaryote results are still scientifically useful

## Proposed scope for the eukaryotic branch

Keep the first version intentionally narrow.

Recommended starting scope:

- **halophilic / saline-adapted algae**
- **high-temperature / thermoacidophilic algae**

Recommended first-pass taxa:

- `Dunaliella salina`
- `Dunaliella parva`
- `Dunaliella tertiolecta`
- `Asteromonas gracilis`
- `Galdieria sulphuraria`
- `Cyanidioschyzon merolae`
- `Cyanidium caldarium`
- `Picochlorum celeri`
- `Picochlorum renovo`

This matches the earlier exploratory algae notes while keeping the restart
manageable.

## Proposed curation rules

### 1. Identity hierarchy

Use a relaxed but explicit identity hierarchy:

1. exact culture-collection strain or isolate
2. named study strain with traceable publication history
3. species-level entry only when phenotype evidence is strong and strain
   mapping is not recoverable

Species-level entries can be acceptable in the exploratory branch, but they
must be marked explicitly as species-level and not treated as equivalent to the
prokaryotic strain-standard core.

### 2. Salinity rules

Do not accept marine habitat words alone as evidence.

Use:

- `saline-adapted alga` for explicit optimum `0.2 <= NaCl equivalent < 0.85 M`
- `halophilic alga` for explicit optimum `>= 0.85 M`

Only accept salinity labels when the paper provides:

- an explicit optimum salinity
- a growth maximum across a salinity gradient
- or a clearly reported optimum range that can be converted defensibly

### 3. Temperature rules

Use a separate algae heat framework:

- `high-temperature alga` for explicit optimum growth `>= 40 C`
- `thermophilic alga` for explicit optimum growth `>= 45 C`
- `thermoacidophilic alga` when high-temperature growth and acidic growth
  context are both explicit

Do not use survival temperatures alone.

### 4. Evidence rules

Only accept a row when:

- the source is citable directly
- the phenotype value is explicit enough to defend
- the identity is stable enough to explain in the thesis

Reject rows based only on:

- habitat inference
- genus reputation
- weak review-language without a usable underlying value

## Proposed technical implementation

The current searchable resource already has fields that support this extension:

- `dataset_branch`
- `dataset_layer`
- `domain_scope`
- `phenotype_focus`

So the cleanest implementation is:

1. keep the current prokaryotic branches unchanged
2. add a new dataset branch such as `eukaryotic_algae_exploratory`
3. set `dataset_layer = exploratory`
4. set `domain_scope = eukaryote`
5. keep these rows excluded from the default ML core unless they pass a later,
   stricter retrieval standard

Recommended new or restarted files:

- `paper_screening_log_eukaryotic_algae_halophilic.tsv`
- `paper_screening_log_eukaryotic_algae_thermophilic.tsv`
- `core_taxon_table_eukaryotic_algae.tsv`
- `sequence_retrieval_tracker_eukaryotic_algae_yes.tsv`

Recommended sequence policy:

1. stabilize metadata first
2. only then attempt proteome retrieval
3. prefer exact strain proteomes
4. fall back to species-level proteomes only with explicit annotation

## What the eukaryotic branch should contribute to the thesis

A realistic contribution is not "a full second core dataset." A better target
is:

- a reproducible eukaryotic curation protocol
- a small, defensible pilot set of accepted taxa
- explicit comparison between prokaryotic and eukaryotic evidence structure
- a statement of what blocks direct unification

That is already thesis-worthy, because it shows where the current framework is
general and where it becomes domain-specific.

## Suggested deliverables

By the end of this branch, aim for:

1. a documented eukaryote-specific curation standard
2. a screened candidate table
3. a small accepted metadata table
4. a retrieval tracker with status labels
5. if possible, a pilot proteome set for a subset of rows
6. integration into the searchable app as a clearly marked exploratory branch

## Suggested success criteria

Call the branch successful if you achieve all of the following:

- at least one clearly documented halophilic/saline-adapted algae subset
- at least one clearly documented high-temperature/thermoacidophilic algae
  subset
- explicit separation between species-level and strain-level entries
- no weakening of the existing prokaryotic core criteria
- searchable integration without double-counting or hidden mixing

## Main risks

The main risks are:

- phenotype evidence is good but sequence linkage is weak
- sequence data exist but are not strain-resolved
- taxonomy/synonym reconciliation becomes too time-consuming
- the branch grows into a second thesis instead of a controlled extension

## Risk control

To keep scope safe:

- cap the first pass at roughly `8-15` accepted eukaryotic rows
- treat this as a pilot branch, not a parity branch with prokaryotes
- keep retrieval as a second phase, not the entry criterion
- stop expansion if the metadata layer is not stable enough

## Recommended wording for supervisor discussion

"I propose adding eukaryotes as a separate exploratory dataset branch rather
than merging them into the current prokaryotic core. The main thesis resource
would remain the strain-aware prokaryotic extremophile dataset, while the
eukaryote branch would test whether the same curation framework can be adapted
to algae with different identity and evidence constraints. This would let us
show framework extensibility without weakening the methodological strictness of
the current core dataset."

## Bottom line

Yes, the thesis can include a eukaryotes branch, but it should be presented as
a **separate exploratory extension** with its own criteria, its own workflow,
and a deliberately limited first-pass scope.
