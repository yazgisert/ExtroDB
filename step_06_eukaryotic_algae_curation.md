# Step 6: Eukaryotic Algae Branch

This branch has been reset and should now be treated as a future separate
dataset, not as part of the current searchable prokaryotic core.

The current thesis core remains:

- prokaryotic thermophiles
- prokaryotic halophiles

Eukaryotic algae should be restarted later from scratch under separate
criteria.

## Why algae need a separate framework

The supervisor-approved direction is to keep algae separate rather than weaken
the prokaryotic core criteria.

That means:

- do not mix algae into the current core dataset
- do not force algae through the same strain strictness as prokaryotes
- do not interpret algae rows as if they were already core-ready

The reason is methodological, not biological doubt:

- algae often have weaker strain linkage in the literature
- culture collection and paper naming are often harder to reconcile
- many relevant algae are thermoacidophilic or hypersaline specialists rather
  than clean matches to the current prokaryotic thresholds

## Proposed algae-specific curation criteria

These criteria are intended for a future algae-specific exploratory branch.
They are separate from the prokaryotic core and should be documented as such.

### 1. Identity requirement

Prefer this identity order:

1. named culture-collection strain or isolate
2. named study strain with traceable publication history
3. species-level entry only if no stable strain mapping exists

For algae, species-level entries may be acceptable in the exploratory branch if
the phenotype evidence is explicit and the limitation is documented.

### 2. Halophilic algae rule

Use `halophilic alga` only when one of these is true:

- an explicit optimum salinity is reported
- a growth maximum is reported across a salinity gradient
- the paper or collection record gives a clear optimum range that can be
  converted to a defensible midpoint

Do not accept habitat words alone such as:

- marine
- brackish
- saline environment
- salt lake isolate

Those are candidate-discovery clues only.

Recommended interpretation for algae:

- `saline-adapted alga`: explicit optimum `0.2 <= NaCl equivalent < 0.85 M`
- `halophilic alga`: explicit optimum `>= 0.85 M`
- keep `moderate` and `extreme` subclasses only when the quantitative evidence
  supports them cleanly

This is stricter than the current prokaryotic halophile threshold because many
ordinary marine eukaryotes already tolerate salinity near seawater conditions.

### 3. Thermophilic algae rule

For algae, use a separate heat-adaptation framework:

- `high-temperature alga`: explicit optimum growth temperature `>= 40 C`
- `thermophilic alga`: explicit optimum growth temperature `>= 45 C`
- `thermoacidophilic alga`: same as above, but with explicit low-pH growth
  context documented in notes

Do not use survival limits or upper tolerance limits by themselves.
Prefer true optimum growth values or reported temperature maxima for active
growth.

### 4. Evidence requirement

Accept an algae row into a future exploratory branch only when:

- the identity is stable enough to explain
- the phenotype value is explicit enough to defend
- the source can be cited directly

Keep it out if the signal is only:

- habitat inference
- genus reputation
- review language without a usable value

### 5. Sequence retrieval rule

Do not start sequence retrieval for algae until the algae metadata table is
stable.

When that time comes, use this retrieval preference:

1. exact strain or isolate proteome/genome proteins
2. culture-collection matched genome
3. species-level proteome only if that fallback is explicit in notes

## Recommended future file strategy

Keep these files, but start them again from empty:

- `paper_screening_log_eukaryotic_algae_halophilic.tsv`
- `paper_screening_log_eukaryotic_algae_thermophilic.tsv`
- `core_taxon_table_eukaryotic_algae.tsv`

## Suggested future starting points

When you restart the algae branch, these are still reasonable directions:

- halophilic algae: `Dunaliella`, `Asteromonas`
- high-temperature or thermoacidophilic algae: `Galdieria`,
  `Cyanidioschyzon`, `Cyanidium`

These are only search directions for the restart, not pre-accepted records.
