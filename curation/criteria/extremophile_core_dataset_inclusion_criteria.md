# Minimum Inclusion Criteria For The First Curated Extremophile TaxID List

This document defines the **minimum bar** for a taxon to enter the first
curated thesis dataset.

The purpose is to keep the initial dataset:
- biologically interpretable
- label-reliable enough for validation tasks
- realistic for a Master's thesis timeline

## Current active scope

At the current stage of the project, these criteria apply to the active
prokaryotic core dataset only:

- thermophile core
- halophile core

Eukaryotic algae are not part of the active core at this time.
They will be restarted later as a separate branch with separate criteria.

## Core principle

Start with a **smaller, phenotype-backed core dataset**.

Candidate-discovery resources such as AFESM can be used to expand coverage, but
they should not automatically define the main training labels.

## Inclusion criteria for the first core list

A taxon should enter the first curated TaxID list only if all of these are met:

### 1. Taxonomy is stable enough

- a stable organism/taxon name is available
- a usable NCBI TaxID can be assigned or reconciled
- taxonomy is specific enough to interpret biologically

Recommended minimum:
- `species` or `genus` preferred
- `family` acceptable only if phenotype support is unusually strong and clearly
  associated

### 2. There is phenotype or strong literature support

At least one of the following should exist:
- explicit phenotype metadata from a trusted source
- a database annotation strongly suggesting extremophile status
- literature support clearly linking the taxon to the target extremophile trait

Examples:
- reported optimal growth temperature
- reported pH optimum / range
- reported salinity optimum / requirement
- repeated literature description as thermophile, halophile, etc.

### 3. The candidate can be mapped to a usable extremophile label

For the first thesis phase, labels should be simple and robust.

Recommended first-pass labels:
- `thermophile`
- `psychrophile`
- `halophile`
- `acidophile`
- `alkaliphile`
- `piezophile`

If the phenotype is too vague to assign one of these labels, keep it out of the
core dataset for now.

For the dedicated cold-side branch, `psychrophile` should be treated as its
own label rather than folded into `thermophile`.

### 4. Evidence is interpretable enough for downstream modeling

You should be able to explain:
- why the taxon was included
- what the label means
- what the evidence source was

If that explanation is weak, the candidate belongs in the exploratory table,
not the main dataset.

## Exclusion criteria for the first core list

Exclude a candidate from the first curated core dataset if any of these apply:

- taxonomy is only broad LCA-level and not biologically specific
- no stable NCBI TaxID can be assigned
- biome information exists, but phenotype support does not
- the taxon is clearly from a non-relevant context for the initial scope
- the organism-label mapping is too ambiguous
- the phenotype is inferred only indirectly from a broad metagenomic source

## Confidence tiers

Use three curation tiers to keep decisions explicit.

### Tier A: include now

All of the following are true:
- stable taxon identity
- usable NCBI TaxID
- clear phenotype support
- clear extremophile label

These taxa go into the first core thesis dataset.

### Tier B: candidate / review

Some useful signal exists, but not enough for direct inclusion.

Examples:
- strong AFESM biome context
- good taxonomy, but missing phenotype values
- literature hints, but not enough for confident labeling

These stay in the candidate table for follow-up enrichment.

### Tier C: exclude for now

Too little support for meaningful inclusion.

Examples:
- taxonomy too unspecific
- no stable taxon mapping
- biome context only
- unclear or conflicting phenotype evidence

## Recommended first-scope strategy

For a Master's thesis, the safest starting point is:

- build a core dataset from phenotype-backed taxa
- choose one or two simpler extremophile categories first
- use broader resources like AFESM as expansion sources

My recommended initial categories:
- `thermophile`
- `halophile`

Why:
- they are easier to motivate biologically
- they are easier to communicate in the proposal
- they are more likely to support interpretable first-pass labels

## What counts as enough evidence in practice

### Strong enough for inclusion

- a taxon with an NCBI TaxID and reported optimum growth temperature
- a taxon described in trusted sources as halophilic with explicit salt
  requirement
- a literature-backed extremophile species with repeated consistent annotation

### Not strong enough for inclusion yet

- a metagenomic cluster from a `Salt` biome with only predicted taxonomy
- a marine biome assignment without clear extremophile phenotype
- an LCA assignment above family level without phenotype evidence

## How to justify this in the proposal

You can say:

“Because downstream validation and model comparison depend on reliable labels, I
plan to start with a phenotype-backed core dataset and treat broader resources
such as AFESM as candidate-discovery and expansion sources rather than direct
ground truth.”

That is the key methodological idea behind these criteria.
