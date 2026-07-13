# Thesis Progress Summary For Supervisor

Last updated: 2026-06-13

## 1. Project aim

The practical goal of the thesis work so far has been to build a curated,
strain-aware extremophile proteome resource that can later be used for:

- comparative biological analysis
- sequence-based machine-learning experiments
- a local searchable web interface for inspection and export

The emphasis has been on defensible phenotype labels rather than collecting
large numbers of sequences without curation.

## 2. What has been completed so far

### 2.1 Operational definitions and criteria

Operational definitions have been formalized and documented for the current
extremophile branches.

Current working thresholds:

- `thermophile`: optimum growth temperature `>= 45 C`
- `hyperthermophile`: optimum growth temperature `>= 80 C`
- `halophile`: NaCl optimum `> 0.2 M`
- `slight halophile`: `0.2 <= NaCl optimum < 0.85 M`
- `moderate halophile`: `0.85 <= NaCl optimum < 3.4 M`
- `extreme halophile`: `NaCl optimum >= 3.4 M`
- `acidophile`: optimum growth `pH <= 3.0`
- `hyperacidophile`: optimum growth `pH <= 1.0`
- `moderate acidophile`: `3.0 < pH <= 5.0`
- `acidotolerant`: `5.0 < pH <= 7.5`

These criteria were not only chosen operationally, but also justified against
review literature and primary strain papers in
[criteria_sources_for_thesis.md](/Users/yazgisert/Documents/Playground/thesis_data_collection/criteria_sources_for_thesis.md).

### 2.2 Curated phenotype branches built so far

The project has progressed beyond the initial thermophile and halophile core.
At the current stage, the curated prokaryotic branches include:

- thermophiles
- halophiles
- acidophiles
- moderate acidophiles
- acidotolerants

For each branch, the workflow now exists from:

1. candidate screening
2. phenotype validation
3. strain-level curation
4. sequence retrieval tracking
5. local protein FASTA integration
6. manifest generation
7. merged FASTA generation
8. searchable index integration

### 2.3 Main curation principle

The curation strategy has consistently prioritized:

- explicit phenotype evidence
- exact or defensible strain mapping
- traceable sequence provenance

Rows were not accepted into the core simply because an organism came from a
relevant habitat. Instead, inclusion required explicit phenotype values or
strong literature/database evidence.

### 2.4 Sequence retrieval and proteome preparation

Accepted rows were linked to sequence sources using a reproducible retrieval
workflow, usually with:

1. exact strain matching where possible
2. GTDB or NCBI assembly/protein accessions
3. explicit documentation of any fallback or ambiguity

Local protein FASTA files were then organized into phenotype-specific raw
folders, and merged FASTA outputs were created for downstream use.

## 3. Current dataset status

### 3.1 Phenotype-specific downloaded proteome sets

Current downloaded branch-level proteome outputs are:

- thermophile: `68 / 69` tracked rows downloaded, `152,071` proteins
- halophile: `46 / 46` tracked rows downloaded, `180,326` proteins
- acidophile: `21 / 21` tracked rows downloaded, `50,761` proteins
- moderate acidophile: `15 / 15` tracked rows downloaded, `55,930` proteins
- acidotolerant: `10 / 10` tracked rows downloaded, `52,214` proteins

Current combined acid-related merged set:

- strict acidophile + moderate acidophile: `36 / 36` rows, `106,691` proteins

### 3.2 Multi-label overlap issue and solution

During expansion, it became clear that some strains fit more than one
extremophile label. For example, some organisms are both thermophilic and
acidophilic.

This created an important methodological issue:

- duplicated proteomes are acceptable for branch-specific curation history
- duplicated proteomes are not acceptable as direct ML input because they can
  inflate counts and introduce leakage

To solve this, the searchable resource was restructured into a hybrid model:

- a deduplicated master proteome index for ML-safe counting
- a separate membership/provenance table that keeps all branch-level label
  assignments

Current overlap count:

- `18` proteomes currently belong to more than one branch membership

### 3.3 Current deduplicated master resource

The current deduplicated master resource contains:

- `145` unique proteomes total
- `142` downloaded proteomes
- `3` proteomes still missing local FASTA files
- `446,497` unique protein sequences in the downloaded master index

This is now the correct ML-facing resource because each proteome is counted
once, even if it has multiple extremophile labels.

## 4. Searchable resource and web app

A local searchable resource and lightweight web app have been built to inspect
the curated dataset.

Current outputs:

- [extremophile_master_index.tsv](/Users/yazgisert/Documents/Playground/thesis_data_collection/searchable_resource/extremophile_master_index.tsv)
- [extremophile_membership_index.tsv](/Users/yazgisert/Documents/Playground/thesis_data_collection/searchable_resource/extremophile_membership_index.tsv)
- [extremophile_master_summary.tsv](/Users/yazgisert/Documents/Playground/thesis_data_collection/searchable_resource/extremophile_master_summary.tsv)

The web app supports:

- filtering by phenotype branch
- filtering by subclass
- filtering by retrieval/download status
- viewing evidence and accession provenance
- opening local FASTA, manifest, and merged files

Importantly, overlapping organisms are now visible under multiple labels in the
app without being counted multiple times in the deduplicated master index.

## 5. Hyperacidophile decision

The optional `hyperacidophile` subclass was explored using the threshold:

- optimum growth `pH <= 1.0`

Conclusion so far:

- biologically the subclass is real
- however, it appears to be genuinely small
- only a very small number of named strains currently satisfy the criterion
  cleanly and are retrieval-ready

Current retrieval-ready hyperacidophile members:

- `Picrophilus oshimae`
- `Picrophilus torridus`

`Ferroplasma thermophilum` was validated biologically as a likely member by
primary literature, but it was not added because an acceptable genome/proteome
source could not currently be retrieved under the project standard.

## 6. Algae branch status

The earlier exploratory algae work has not been merged into the active
prokaryotic resource.

Current position:

- algae remain outside the active searchable core
- the algae branch is effectively a future restart branch
- it will require separate curation logic and source hierarchy

This separation was important because the algae workflow needs a different
identity and phenotype-evidence standard than the prokaryotic branches.

## 7. Main methodological achievements so far

The most important completed methodological advances are:

1. A reproducible curation framework for extremophile strain selection.
2. Operational phenotype thresholds that are explicit and literature-backed.
3. A sequence retrieval workflow tied to curated strain metadata.
4. Branch-specific raw FASTA, manifest, and merged FASTA organization.
5. A searchable local resource for auditing and export.
6. A proteome-deduplicated master index suitable for downstream ML work.

## 8. Current limitations

The main limitations at the current stage are:

- a few thermophile rows are still missing local FASTA files
- the acidophile branches are smaller than the halophile branch
- the hyperacidophile subclass is biologically interesting but numerically
  sparse
- the algae branch has not yet been rebuilt under its separate framework

## 9. Recommended next steps

The most logical next steps are:

1. resolve the remaining missing thermophile proteomes
2. decide whether to expand thermophile coverage further for better balance
3. use the deduplicated master proteome index as the ML-facing resource
4. define the final export format for model-ready multi-label training data
5. restart the algae branch later under a clearly separated workflow

## 10. Short thesis-ready summary sentence

At this stage, the thesis has produced a curated, strain-aware extremophile
proteome resource spanning thermophiles, halophiles, acidophiles, moderate
acidophiles, and acidotolerants, together with a deduplicated searchable master
index that preserves phenotype overlap while preventing proteome double-counting
for downstream machine-learning applications.
