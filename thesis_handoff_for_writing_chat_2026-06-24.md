# Thesis Handoff For Writing/Structuring Chat

Last updated: 2026-06-24

## 1. What this document is for

This document is meant to be uploaded into another chat so that the next model
has a clear, detailed understanding of:

- what the thesis is about
- what has already been built
- how the dataset and codebase are organized
- what methodological decisions have already been made
- what was updated most recently
- what still remains to be written, analyzed, or modeled

The goal is to make it easy for another chat to help structure the thesis
writing without needing to rediscover the whole project from scratch.

## 2. Thesis title and core direction

Current working thesis title:

- `Specializing Protein Language Models for Protein Engineering using Extremophiles`

In plain language, the thesis is about building a biologically defensible,
curated extremophile protein/proteome dataset and then using that dataset to
test whether protein language models (PLMs) benefit from specialization on
extremophile data.

This means the thesis is not only a machine-learning thesis and not only a
biology curation thesis. It sits at the intersection of:

- biological data curation
- sequence dataset construction
- representation learning on proteins
- evaluation of domain specialization for PLMs
- future support for protein engineering workflows

## 3. The core thesis problem

The central problem is that extremophile proteins are biologically interesting
for protein engineering, but the available public data are fragmented,
inconsistently labeled, and often not clean enough for machine-learning use.

More specifically:

- many public resources do not provide thesis-grade phenotype labels
- phenotype evidence is often incomplete or only indirectly supported
- strain identity can be ambiguous
- sequence sources may be available, but not clearly tied to the curated strain
- overlapping extremophile categories create leakage risks for downstream ML

So the thesis is trying to solve a foundation problem:

- first create a curated, phenotype-backed, sequence-linked extremophile
  dataset
- then test whether this dataset can support meaningful PLM specialization

## 4. The thesis aim

The current practical aim of the thesis can be stated like this:

Build a curated extremophile-focused dataset with explicit phenotype evidence,
stable strain-aware sequence mapping, and reproducible provenance, then use it
as the basis for evaluating whether protein language models can be specialized
usefully on extremophile data.

A slightly more detailed version:

1. Collect and curate extremophile organisms using explicit phenotype rules.
2. Link accepted rows to defensible proteome or genome-protein sequence sets.
3. Organize the resulting data into branch-specific and master resources.
4. Create a searchable local resource for inspection, auditing, and export.
5. Use the cleaned dataset as the ML-facing basis for later PLM analysis and
   fine-tuning.

## 5. The core hypothesis

The central hypothesis is:

Protein language models specialized on a curated extremophile dataset may learn
more domain-relevant sequence representations than fully generalist PLMs.

Important nuance:

- this is treated as a testable hypothesis, not an assumption
- a negative result would still be scientifically useful
- even if specialization does not improve performance, the thesis still
  produces value through the curated dataset, explicit methodology, and
  evaluation-ready resource layer

## 6. Current research-question framing

The thesis has been framed around questions like these:

1. Can a phenotype-backed extremophile dataset be assembled in a reliable and
   interpretable way for PLM specialization?
2. Does the resulting dataset contain enough meaningful signal for quality
   checks and downstream predictive tasks?
3. Does PLM specialization on extremophile data improve useful representations
   or downstream behavior compared with generalist PLMs?
4. Can exploratory large-scale resources be used for candidate discovery
   without compromising label quality?
5. Do current generalist PLMs appear to model extremophile proteins less well
   than more typical proteins, for example in pseudoperplexity-style analyses?

## 7. What has been accomplished so far at a high level

The major completed work so far is not yet the final PLM modeling stage. The
main completed body of work is the full dataset-construction and resource layer
that the later modeling depends on.

In practice, the completed achievements so far are:

1. Operational phenotype definitions were formalized.
2. Candidate-screening workflows were created for multiple extremophile
   branches.
3. Strain-aware curation tables were built.
4. Sequence retrieval trackers were built and populated.
5. Local raw proteome FASTA collections were organized by branch.
6. Manifest tables and merged FASTA outputs were generated.
7. A searchable master dataset was built.
8. A local web app was built to inspect and export the resource.
9. Overlap handling was improved so the ML-facing dataset is deduplicated at
   the proteome level.

This is important for thesis writing because the work is already substantial
even before model fine-tuning begins. The thesis already includes a real
methods contribution in dataset design, curation logic, provenance handling,
and resource building.

## 8. The biological scope of the active dataset

The current active searchable dataset is focused on phenotype-backed
prokaryotic extremophile entries.

Active branches currently represented in the live resource are:

- `psychrophile_core`
- `thermophile_core`
- `halophile_core`
- `acidophile_expansion`
- `moderate_acidophile_expansion`
- `acidotolerant_expansion`
- `alkaliphile_core`

Important note:

- some older summary documents in the repository predate the current full live
  state and may underreport newer branches such as `alkaliphile_core`
- the live `searchable_resource` files are the most current representation of
  the active dataset state

## 9. Current phenotype definitions

The main operational rules currently in use are:

### 9.1 Temperature-related labels

- `psychrophile`: `ogt_celsius <= 15`
- `thermophile`: `ogt_celsius >= 45`
- `hyperthermophile`: `ogt_celsius >= 80`

Excluded:

- thermotolerant organisms with `40 <= ogt_celsius < 45`
- psychrotolerant/psychrotrophic organisms if `ogt_celsius > 15`

### 9.2 Halophile labels

- `halophile`: `nacl_optimum_M > 0.2`
- `slight halophile`: `0.2 <= nacl_optimum_M < 0.85`
- `moderate halophile`: `0.85 <= nacl_optimum_M < 3.4`
- `extreme halophile`: `nacl_optimum_M >= 3.4`

Excluded:

- halotolerant organisms without an optimum above `0.2 M`

### 9.3 Acid-related labels

- `acidophile`: `ph_optimum <= 3.0`
- `hyperacidophile`: `ph_optimum <= 1.0`
- `moderate acidophile`: `3.0 < ph_optimum <= 5.0`
- `acidotolerant`: `5.0 < ph_optimum <= 7.5`

### 9.4 Alkaliphile labels

- `alkaliphile`: `ph_optimum >= 9.0`
- `extreme alkaliphile`: `ph_optimum >= 10.0`
- `alkalitolerant`: `7.5 < ph_optimum < 9.0` with documented growth at
  `pH >= 9.0`

## 10. Curation philosophy

The most important methodological principle is that the dataset is not built by
simply collecting organisms from extreme habitats.

Instead, inclusion aims to require:

- explicit phenotype evidence
- defensible identity mapping
- stable sequence provenance
- traceable justification in notes or evidence fields

This means:

- habitat implication alone is not enough
- species-only evidence is usually weaker than strain-resolved evidence
- exact strain mapping is preferred whenever possible
- every accepted row should be citable and defendable in a thesis context

## 11. Evidence and inclusion logic

The dataset uses evidence categories such as:

- `experimental`
- `curated_database`
- `review`
- `genus_inference`
- `biome_inference`

The inclusion logic has historically allowed values like:

- `include_in_core = yes`
- `include_in_core = review`
- `include_in_core = no`

However, the current live searchable dataset has been cleaned so that there are
no remaining review-status proteome rows in the active resource. The recent UI
and data cleanup specifically removed the need to expose a "proteome status"
review filter in the web app.

## 12. Data workflow in plain language

The current workflow can be described simply as:

1. Define phenotype criteria.
2. Screen candidate taxa from literature or databases.
3. Validate phenotype evidence.
4. Curate accepted organism/strain rows into branch-specific core tables.
5. Build sequence retrieval trackers for accepted rows.
6. Link rows to GTDB or NCBI sequence sources.
7. Download or place local protein FASTA files.
8. Build manifests summarizing available local FASTA coverage.
9. Build merged FASTA files per phenotype branch.
10. Build a searchable master resource for auditing and export.
11. Use the deduplicated master resource for ML-facing downstream work.

## 13. Main data objects in the project

The key file types are:

### 13.1 Screening logs

Purpose:

- first-pass evidence collection
- candidate tracking
- phenotype-value extraction
- preliminary accept/reject/review decisions

Examples:

- `paper_screening_log_psychrophile.tsv`
- `paper_screening_log_tempura_bacdive.tsv`
- `paper_screening_log_halophile_bacdive.tsv`
- `paper_screening_log_acidophile_bacdive.tsv`

### 13.2 Core taxon tables

Purpose:

- curated organism-level rows
- phenotype flags and optimum values
- evidence metadata
- inclusion decisions

Examples:

- `core_taxon_table_psychrophile.tsv`
- `core_taxon_table_tempura_bacdive_yes.tsv`
- `core_taxon_table_halophile_bacdive.tsv`
- `core_taxon_table_acidophile_bacdive.tsv`
- `core_taxon_table_moderate_acidophile_bacdive.tsv`
- `core_taxon_table_acidotolerant_bacdive.tsv`

### 13.3 Sequence retrieval trackers

Purpose:

- map curated rows to real sequence sources
- record GTDB/NCBI accessions
- document retrieval status and fallback logic

Examples:

- `sequence_retrieval_tracker_psychrophile_yes.tsv`
- `sequence_retrieval_tracker_tempura_bacdive_yes.tsv`
- `sequence_retrieval_tracker_halophile_bacdive_yes.tsv`
- `sequence_retrieval_tracker_acidophile_bacdive_yes.tsv`

### 13.4 Raw protein FASTA folders

Purpose:

- store the branch-specific local proteome files

Examples:

- `psychrophile_raw_protein_fastas/`
- `thermophile_raw_protein_fastas/`
- `halophile_raw_protein_fastas/`
- `acidophile_raw_protein_fastas/`
- `moderate_acidophile_raw_protein_fastas/`
- `acidotolerant_raw_protein_fastas/`
- `alkaliphile_raw_protein_fastas/`

### 13.5 Manifest tables

Purpose:

- summarize which tracker rows currently have local FASTA files
- record protein counts
- connect trackers to local sequence assets

### 13.6 Merged FASTA files

Purpose:

- create branch-level concatenated proteome sets for downstream analysis

### 13.7 Searchable resource tables

Purpose:

- provide the app-facing and ML-facing structured dataset layer

Key files:

- `searchable_resource/extremophile_master_index.tsv`
- `searchable_resource/extremophile_membership_index.tsv`
- `searchable_resource/extremophile_master_summary.tsv`

## 14. Why the searchable resource matters

The searchable resource is an important thesis contribution because it turns
many separate curation outputs into a usable local resource.

It helps with:

- auditing phenotype assignments
- checking evidence provenance
- inspecting sequence source links
- opening local FASTA and manifest files
- exporting filtered subsets
- making the dataset easier to describe in the thesis methods/results sections

Conceptually, it also marks the transition from "many branch-specific working
tables" to "one coherent resource layer."

## 15. The overlap problem and why it mattered

A major methodological issue appeared once multiple phenotype branches became
active:

- some organisms fit more than one extremophile category
- for example, an organism can be both thermophilic and acidophilic

This matters because:

- overlap is biologically real
- but naive counting would duplicate proteomes
- duplicated proteomes can distort statistics
- duplicated proteomes are especially problematic for downstream ML because
  they can create data leakage

The solution used here was a hybrid design:

- a deduplicated master proteome index for ML-facing work
- a separate membership/provenance layer that preserves multi-branch
  assignments

This is one of the strongest methodological points in the thesis, because it
shows awareness of both biological reality and ML evaluation hygiene.

## 16. Current live dataset state as of 2026-06-24

These counts reflect the current live `extremophile_master_index.tsv` after the
most recent cleanup.

### 16.1 Master index totals

- unique master rows: `213`
- summed proteins across the master index: `713,417`
- rows with more than one branch membership: `18`

Important interpretation:

- the master index is deduplicated at the proteome-entry level
- branch memberships can still overlap
- therefore branch counts do not sum cleanly to the master total

### 16.2 Current branch-level totals in the live resource

- `acidophile_expansion`: `21` rows, `50,761` proteins
- `acidotolerant_expansion`: `10` rows, `52,214` proteins
- `alkaliphile_core`: `43` rows, `149,056` proteins
- `halophile_core`: `46` rows, `180,326` proteins
- `moderate_acidophile_expansion`: `15` rows, `55,930` proteins
- `psychrophile_core`: `28` rows, `117,864` proteins
- `thermophile_core`: `68` rows, `152,071` proteins

### 16.3 Current phenotype-label totals across the live master index

These are label totals, not branch totals:

- `psychrophile = yes`: `28`
- `thermophile = yes`: `71`
- `hyperthermophile = yes`: `32`
- `halophile = yes`: `58`
- `acidophile = yes`: `36`
- `alkaliphile = yes`: `29`

Important distinction:

- branch totals count membership in named branches such as `halophile_core`
- label totals count the presence of a phenotype flag
- those numbers can differ because some rows have multiple labels or belong to
  differently named branch structures

## 17. What was done recently in the web app and dataset

The most recent work in the current chat session included both UI cleanup and
dataset cleanup.

### 17.1 Detail-panel UI fix

Problem:

- long detail values in the right-hand detail section were visually spilling
  outside the main panel

What was changed:

- the detail card and its items were given safer width and wrapping behavior
- long strings such as notes and local FASTA paths now wrap inside the panel

### 17.2 Filter simplification

Problem:

- the filter section still exposed controls for `Dataset Branch` and
  `Proteome Status`
- these were no longer aligned with the desired current use case

What was changed:

- the `Dataset Branch` filter was removed from the UI
- the `Proteome Status` filter was removed from the UI
- the associated JavaScript filter logic was removed as well

### 17.3 Review-status check

Question checked:

- whether any `review` proteome statuses were still present in the live data

Result:

- there were no `review` proteome-status rows in the live searchable dataset

### 17.4 Removal of the remaining missing entries

Three remaining missing rows were identified and removed from the live
searchable resource:

- `Actibacter sediminis` (`JC2129`) from `halophile_core`
- `Actinoallomurus coprocola` (`TT04-09`) from `halophile_core`
- `Methanobacterium thermaggregans` (`DSM 3266`) from `thermophile_core`

Files updated:

- `searchable_resource/extremophile_master_index.tsv`
- `searchable_resource/extremophile_membership_index.tsv`
- `searchable_resource/extremophile_master_summary.tsv`

Effect:

- `halophile_core` now contains `46` entries total
- `thermophile_core` now contains `68` entries total
- the active live searchable dataset now no longer contains those unresolved
  missing rows

### 17.5 Web app recoloring

The app theme was also updated visually:

- the top hero/header region was changed to a purple gradient using
  `#9D19D1` and `#B774D7`
- the previous green accent color was changed to the purple palette throughout
  the main app accents

This is only a presentation/UI change and does not affect the dataset logic.

## 18. What the thesis has not completed yet

Although a lot has been built, the thesis is not finished.

The largest unfinished parts are on the ML-analysis and final-writing side.

Still incomplete or not yet finalized:

- final full thesis manuscript
- final polished introduction/background chapters
- final methods write-up in thesis form
- final descriptive results section with tables/figures
- baseline predictive evaluation and/or dataset-quality experiments
- PLM pseudoperplexity or related representation analyses
- PLM fine-tuning and comparison of specialized vs generalist models
- final discussion and limitations chapter

So the thesis currently has a strong data foundation, but the modeling and
final narrative still need to be completed.

## 19. What the thesis already contributes even before modeling

If another chat helps structure the thesis, it should not frame the work as if
the real contribution only begins with fine-tuning.

The thesis already contributes:

1. A curation framework for phenotype-backed extremophile data.
2. Explicit operational definitions for multiple extremophile classes.
3. A sequence-linking pipeline tied to curated metadata.
4. Branch-specific proteome sets and merged FASTA outputs.
5. A deduplicated master resource designed with ML leakage concerns in mind.
6. A local searchable resource/web app for audit and export.

That means the thesis can be written as a methods-heavy contribution even if
the final modeling section remains smaller than originally envisioned.

## 20. Likely thesis narrative

A strong writing narrative would probably look something like this:

### 20.1 Problem

Extremophile proteins are interesting for protein engineering, but there is no
clean, phenotype-backed, ML-ready extremophile benchmark/resource.

### 20.2 Gap

Public data are fragmented, labels are inconsistent, and sequence provenance is
often not explicit enough for reliable ML work.

### 20.3 Thesis response

Build a curated extremophile resource with explicit phenotype rules, traceable
evidence, and defensible sequence linkage.

### 20.4 Dataset/method contribution

Create branch-level curated sets, manifests, merged FASTAs, and a deduplicated
master index that preserves biological overlap without introducing ML leakage.

### 20.5 Resource contribution

Build a searchable interface to inspect and export the dataset.

### 20.6 ML contribution

Use the curated dataset as the basis for testing whether extremophile-focused
PLM specialization is useful.

## 21. Suggested thesis chapter logic

This is not mandatory, but it is a sensible structure:

1. Introduction
2. Background
3. Data sources and curation framework
4. Dataset construction and sequence-linking workflow
5. Searchable resource and final dataset organization
6. Dataset characteristics and descriptive results
7. PLM evaluation and/or specialization experiments
8. Discussion
9. Conclusion and future work

More detailed interpretation:

### Chapter 1: Introduction

Should explain:

- why extremophiles matter
- why PLM specialization is relevant
- why curated data are the bottleneck
- what the thesis aims to do

### Chapter 2: Background

Should explain:

- protein language models
- extremophile biology
- dataset quality and benchmark issues
- why phenotype-backed data matter
- why overlap and leakage matter for ML

### Chapter 3: Data Sources and Curation Framework

Should explain:

- literature and database sources
- operational definitions
- evidence types
- inclusion/exclusion rules
- strain identity logic

### Chapter 4: Dataset Construction and Sequence Linking

Should explain:

- screening logs
- core tables
- retrieval trackers
- GTDB/NCBI linking
- local FASTA generation
- manifests and merged FASTAs

### Chapter 5: Searchable Resource Layer

Should explain:

- why the searchable layer was built
- master index vs membership index
- overlap handling
- web app functionality
- export and audit use cases

### Chapter 6: Dataset Characteristics

Should explain:

- branch sizes
- phenotype class composition
- overlap statistics
- protein-count totals
- what the live dataset currently contains

### Chapter 7: Modeling / Evaluation

Should explain:

- what PLM analyses were run
- what task framing was used
- whether specialization helped
- what counts as success or informative failure

### Chapter 8: Discussion

Should explain:

- what worked
- what remains limited
- what the dataset is strong at
- where sparsity, ambiguity, or imbalance remain
- what the implications are for future extremophile ML

### Chapter 9: Conclusion

Should explain:

- what the thesis achieved overall
- why the curated resource matters
- how it supports future protein engineering/modeling work

## 22. Important limitations to mention honestly

The thesis should state limitations clearly.

Key limitations include:

- some phenotype categories are still relatively small
- branch sizes are imbalanced
- phenotype evidence quality is not uniform across all taxa
- some earlier planning documents lag behind the current live dataset state
- the final ML specialization/evaluation stage may still be narrower than the
  original broad vision
- future engineering integration is out of scope as a direct thesis
  deliverable

These are not fatal problems. They are normal thesis-scope realities and can
be presented honestly.

## 23. Recommended framing for another chat

If you upload this to another chat and ask it to help structure your thesis
writing, the next chat should understand that:

- the thesis is already well advanced on the data/method side
- the biggest writing task is to turn the curation and resource-building work
  into a coherent thesis narrative
- the master story is "curated extremophile resource first, PLM
  specialization second"
- the thesis should not be framed as if the dataset were only a preliminary
  step

## 24. Short summary that could be pasted into another chat

I am writing a Master's thesis titled "Specializing Protein Language Models for
Protein Engineering using Extremophiles." The work so far has focused on
building a curated, phenotype-backed, strain-aware extremophile proteome
resource rather than jumping directly into modeling. I defined operational
criteria for psychrophiles, thermophiles, halophiles, acidophiles, moderate
acidophiles, acidotolerants, and alkaliphiles; built screening logs, curated
core tables, sequence retrieval trackers, raw FASTA collections, manifests,
merged FASTAs, and a searchable master resource with a local web app. A major
methodological issue was overlap between extremophile categories, so the
dataset was redesigned into a deduplicated master proteome index plus a
separate membership/provenance layer to avoid ML leakage while preserving
multi-label biology. As of 2026-06-24, the live master index contains 213
unique rows and 713,417 proteins, with active branch totals of 28
psychrophiles, 68 thermophile-core entries, 46 halophile-core entries, 21
acidophiles, 15 moderate acidophiles, 10 acidotolerants, and 43
alkaliphile-core entries, plus 18 overlapping multi-branch members. In the
most recent cleanup, I fixed UI overflow in the app, removed unnecessary
branch/status filters, removed 3 remaining missing entries from the live
dataset, and recolored the app header/accent theme. What I need now is help
structuring the thesis narrative and chapters around the idea that the thesis's
first major contribution is the curated extremophile dataset/resource and the
second major contribution is the later evaluation of PLM specialization on that
resource.
