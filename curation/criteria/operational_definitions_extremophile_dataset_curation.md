# Operational Definitions For Extremophile Dataset Curation

This document records the active curation rules for the thesis dataset at the
current stage of the project.

## Current active scope

The active searchable dataset currently includes the phenotype-backed
prokaryotic core and its current pH-related expansions:

- psychrophile core
- thermophile core
- halophile core
- acidophile expansion
- moderate acidophile expansion
- acidotolerant expansion

The psychrophile branch now has an expanded curated core and a populated
sequence-retrieval tracker for accepted rows.

Eukaryotic algae are not part of the active core at this time.
They will be restarted later as a separate branch with separate criteria.

## Current unit of curation

Each row should represent a curated organism entry, ideally at the strain
level.

For the current prokaryotic core:

- strain-level mapping is preferred
- species-level-only evidence is usually not enough for `include_in_core = yes`
- accepted core rows should support sequence retrieval to a matched proteome or
  genome-protein set when possible

## Current prokaryotic phenotype definitions

### Temperature-related labels

The hot-side thermophile labels are active in the current core.
The psychrophile rule below is active in the separate cold-temperature branch.

- `thermophile`: `ogt_celsius >= 45`
- `hyperthermophile`: `ogt_celsius >= 80`
- `psychrophile`: `ogt_celsius <= 15`

Exclude thermotolerant organisms with:

- `40 <= ogt_celsius < 45`

Exclude psychrotolerant or psychrotrophic organisms from a strict
psychrophile branch when:

- low-temperature growth is reported
- but `ogt_celsius > 15`

Treat psychrophiles as a separate cold-temperature branch or as a sibling
branch within a broader `temperature-related` umbrella. Do not store them
inside a branch that remains explicitly named `thermophile`.

### Halophile labels

- `halophile`: `nacl_optimum_M > 0.2`
- `slight halophile`: `0.2 <= nacl_optimum_M < 0.85`
- `moderate halophile`: `0.85 <= nacl_optimum_M < 3.4`
- `extreme halophile`: `nacl_optimum_M >= 3.4`

Exclude halotolerant organisms when the optimum does not exceed `0.2 M`.

### Acid-related labels

- `acidophile`: `pH_optimum <= 3.0`
- `hyperacidophile`: `pH_optimum <= 1.0`
- `moderate acidophile`: `3.0 < pH_optimum <= 5.0`
- `acidotolerant`: `5.0 < pH_optimum <= 7.5`

Exclude acid-tolerant organisms from the strict acidophile branch when the
optimum is above `pH 3.0`.

### Alkaliphile labels

- `alkaliphile`: `pH_optimum >= 9.0`
- `extreme alkaliphile`: `pH_optimum >= 10.0`
- `alkalitolerant`: `7.5 < pH_optimum < 9.0`, with documented growth at
  `pH >= 9.0`

Exclude organisms with only neutral or weakly alkaline optima from the
alkaliphile core, even if they were isolated from alkaline habitats.

## Current evidence definitions

Use these evidence types:

- `experimental`
- `curated_database`
- `review`
- `genus_inference`
- `biome_inference`

## Current inclusion decisions

### `include_in_core = yes`

Use `yes` only when:

- identity is stable enough to defend
- the phenotype value or label is explicit enough to defend
- the threshold is satisfied
- the evidence source is directly citable
- strain ambiguity is low enough for a thesis-grade core row

### `include_in_core = review`

Use `review` when:

- the candidate is biologically promising
- but the value is borderline, incomplete, or conflicting
- strain mapping is unclear
- only species-level evidence is available

### `include_in_core = no`

Use `no` when:

- the threshold is not satisfied
- the source only implies habitat
- the label cannot be defended from the available evidence

## Current sequence-linking rule

For accepted core rows, use this retrieval order:

1. GTDB match to the curated strain
2. NCBI Assembly match when GTDB is not available

Preferred retrieval unit:

1. exact strain proteome
2. exact strain genome proteins
3. species-level fallback only if this is explicitly documented

## Searchable resource rule

The searchable master index should currently include only the active
prokaryotic core.

It should remain filterable by:

- phenotype category
- subclass
- evidence type
- inclusion decision
- retrieval status
- sequence source

## Eukaryotic algae status

The earlier algae records have been cleared.
Algae are now treated as a future separate branch that will restart from empty
tables.

The reason is not that algae are irrelevant.
The reason is that they need a different curation standard than the current
prokaryotic core.

## Proposed future algae framework

When algae are restarted, use a separate exploratory framework.

This branch should remain separate from the prokaryotic core at every stage:

- separate screening tables
- separate curated metadata table
- separate sequence-tracking tables
- no automatic inclusion in the current searchable prokaryotic resource

### Algae candidate-discovery sources

Use source types in this order:

1. `AlgaeBase` for candidate discovery, accepted names, and synonym cleanup
2. culture-collection records such as `UTEX`, `CCAP`, and `NCMA` for
   strain-level identity anchoring
3. primary literature for explicit phenotype values
4. review articles only as supporting context, not as the sole basis for
   thesis-grade acceptance when a primary source should exist

`AlgaeBase` should be treated mainly as a discovery and taxonomic
normalization layer.
It can help identify candidate taxa, accepted names, and known synonyms, but
it should not by itself define the phenotype label.

Culture-collection records should be treated mainly as an identity layer.
They are especially useful for:

- matching strain names across papers
- stabilizing isolate or accession names
- confirming whether a named strain is a real collection-backed organism
- supporting later sequence retrieval at the exact strain or isolate level

### Algae identity rule

For the restarted algae branch, use this identity preference order:

1. named culture-collection strain or isolate
2. named study strain with traceable publication history
3. species-level entry only when no stable strain mapping is available

Species-level entries may still be acceptable in the exploratory algae branch,
but only when:

- the phenotype evidence is explicit
- the identity limitation is documented in `notes`
- the row is not presented as strain-resolved when it is not

### Algae phenotype-evidence rule

Habitat words alone are not enough.

Do not accept a row on the basis of:

- `marine`
- `brackish`
- `saline environment`
- `hot spring`
- `acid spring`

Those are candidate-discovery clues only.

Prefer explicit physiological evidence such as:

- optimum growth salinity
- optimum growth temperature
- experimentally supported growth maximum across a salinity gradient
- experimentally supported active growth at high temperature
- explicit acidic growth context when using `thermoacidophilic alga`

### Recommended future algae labels

- `saline-adapted alga`: explicit salinity optimum `0.2 <= x < 0.85 M`
- `halophilic alga`: explicit salinity optimum `>= 0.85 M`
- `high-temperature alga`: explicit growth optimum `>= 40 C`
- `thermophilic alga`: explicit growth optimum `>= 45 C`
- `thermoacidophilic alga`: above temperature rule plus explicit acidic growth
  context
- `psychrophilic alga`: explicit optimal growth temperature `<= 15 C`, or
  explicit culture-maintenance temperature used as a documented proxy when no
  formal OGT experiment is available and that caveat is stated in `notes`
- `acidophilic alga`: explicit optimal growth pH `<= 3.0` without a
  thermophilic component; if both acidity and high-temperature growth are
  present, prefer `thermoacidophilic alga`

Keep `moderate` and `extreme` halophile subclasses for algae only when the
quantitative salinity evidence supports those boundaries cleanly.

### Algae inclusion rule

For the restarted exploratory algae branch:

- `include_in_core = yes` may be used only within the algae branch itself, not
  the active prokaryotic core
- accepted algae rows must have explicit phenotype evidence and defensible
  identity mapping
- `review` should be used when phenotype evidence is promising but identity,
  thresholds, or source quality remain incomplete
- `no` should be used when the evidence is only habitat-based, contradictory,
  or below threshold

### Algae sequence-linking rule

Do not begin sequence retrieval until the algae metadata table is stable.

When sequence retrieval begins, use this preference order:

1. exact strain or isolate proteome/genome proteins
2. culture-collection matched genome or proteome
3. species-level fallback only if the fallback is explicit in `notes`

This algae workflow is exploratory and literature-informed, but it must still
be explicit and reproducible.
It should be reported as a separate dataset layer from the active
prokaryotic core.
