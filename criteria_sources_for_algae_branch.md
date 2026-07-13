# Sources For Criteria And Labels In The Eukaryotic Algae Branch

This note is the algae-branch counterpart to the prokaryotic criteria file.

Its purpose is to separate:

1. source-backed biological statements
2. project-operational curation rules

That distinction matters even more for algae than for the prokaryotic core,
because some of the most useful algae thresholds are not universal community
standards. They are explicit thesis choices made to keep the branch
reproducible and biologically defensible.

## 0. Current branch status

The algae branch is a separate exploratory branch.

It is **not** part of the active prokaryotic searchable core.

That separation is documented in:

- [step_06_eukaryotic_algae_curation.md](/Users/yazgisert/Documents/Playground/thesis_data_collection/step_06_eukaryotic_algae_curation.md)
- [operational_definitions_extremophile_dataset_curation.md](/Users/yazgisert/Documents/Playground/thesis_data_collection/operational_definitions_extremophile_dataset_curation.md)
- [algae_branch_screening_workflow_and_templates.md](/Users/yazgisert/Documents/Playground/thesis_data_collection/algae_branch_screening_workflow_and_templates.md)

The current intended workflow is:

1. taxonomic cleanup
2. identity reconciliation
3. phenotype evidence screening
4. accepted metadata table
5. sequence retrieval only after metadata is stable

## 1. Source hierarchy for the algae branch

### 1.1 Taxonomy and synonym cleanup

Use `AlgaeBase` as the first discovery and accepted-name layer.

Use it for:

- accepted names
- synonym cleanup
- taxonomic direction

Do **not** use it by itself as phenotype proof.

### 1.2 Identity anchoring

Use culture-collection resources as the identity layer.

Verified official resources:

- [CCAP](https://www.ccap.ac.uk/)
- [NCMA](https://ncma.bigelow.org/)

Why these matter:

- CCAP states that its catalogue is a searchable database containing strain
  data, including origin, maintenance information, references, and sequencing
  data.
- NCMA provides public strain-oriented browsing, including algae categories,
  hypersaline categories, tree-of-life navigation, and strain/publication
  discovery pages.

These are useful for thesis-grade identity support, but they should not be
treated as phenotype ground truth unless the phenotype is explicitly stated in
the collection record itself.

### 1.3 Phenotype evidence

Use primary literature first.

Verified example source pages currently used in this branch:

- [Osmoregulation in the Halotolerant Alga Asteromonas gracilis](https://pubmed.ncbi.nlm.nih.gov/16661724/)
- [Bioenergy application of Dunaliella salina SA 134 grown at various salinity levels for lipid production](https://pubmed.ncbi.nlm.nih.gov/28808229/)
- [Long-term hyposaline and hypersaline stresses produce distinct antioxidant responses in the marine alga Dunaliella tertiolecta](https://pubmed.ncbi.nlm.nih.gov/14610888/)
- [Growth under Different Trophic Regimes and Synchronization of the Red Microalga Galdieria sulphuraria](https://pubmed.ncbi.nlm.nih.gov/34202768/)
- [The carbon-concentrating mechanism of the extremophilic red microalga Cyanidioschyzon merolae](https://pmc.ncbi.nlm.nih.gov/articles/PMC10154280/)
- [Temperature-Induced Remodeling of the Photosynthetic Machinery Tunes Photosynthesis in the Thermophilic Alga Cyanidioschyzon merolae](https://pmc.ncbi.nlm.nih.gov/articles/PMC5411153/)

Use review articles only as supporting context.

## 2. Source-backed biology relevant for the branch

### 2.1 Strong salinity evidence example: Asteromonas gracilis

The PubMed abstract for `Asteromonas gracilis` explicitly states:

- growth occurs from `0.5 M NaCl` to saturation
- the specific growth rate is maximal between `0.5 and 2.5 M NaCl`

Source:

- [PubMed 16661724](https://pubmed.ncbi.nlm.nih.gov/16661724/)

This is exactly the type of explicit salinity-growth evidence that should be
accepted for the algae branch.

### 2.2 Strong strain-linked salinity evidence example: Dunaliella salina

The verified PubMed page confirms the exact study focus:

- `Dunaliella salina SA 134`
- growth under various salinity levels

Source:

- [PubMed 28808229](https://pubmed.ncbi.nlm.nih.gov/28808229/)

This is a good example of a source that is useful because it is already
strain-specific and salinity-focused.

### 2.3 Thermoacidophilic red algae are real and relevant

The verified PubMed page for `Galdieria sulphuraria` explicitly describes it
as an `extremophilic` red microalga.

Source:

- [PubMed 34202768](https://pubmed.ncbi.nlm.nih.gov/34202768/)

The verified PMC page for `Cyanidioschyzon merolae` explicitly describes it as
an `extremophilic red microalga` that grows in `low-pH, high-temperature`
environments.

Source:

- [PMC10154280](https://pmc.ncbi.nlm.nih.gov/articles/PMC10154280/)

The verified PMC temperature-acclimation paper provides further direct context
for `C. merolae` as a thermophilic alga.

Source:

- [PMC5411153](https://pmc.ncbi.nlm.nih.gov/articles/PMC5411153/)

Together, these sources support the thesis claim that a eukaryotic
thermoacidophilic algae branch is biologically justified.

## 3. What is not yet a universal algae standard

This branch should **not** pretend that all of its thresholds are universal
literature definitions in the way some prokaryotic extremophile definitions
are.

Instead, the clean thesis position is:

- some source pages directly justify the candidate phenotype
- the exact branch thresholds are then chosen operationally for reproducible
  curation

This is especially important for:

- saline versus truly halophilic marine algae
- high-temperature versus thermophilic algae
- strain-level versus species-level acceptability

## 4. Identity rule for algae

### Active project rule

Use this identity preference order:

1. named culture-collection strain or isolate
2. named study strain with traceable publication history
3. species-level entry only when no stable strain mapping exists

This is a **project-operational rule**.

It is documented in:

- [operational_definitions_extremophile_dataset_curation.md](/Users/yazgisert/Documents/Playground/thesis_data_collection/operational_definitions_extremophile_dataset_curation.md)
- [step_06_eukaryotic_algae_curation.md](/Users/yazgisert/Documents/Playground/thesis_data_collection/step_06_eukaryotic_algae_curation.md)

### Thesis-ready interpretation

For algae, species-level entries may still be acceptable in the exploratory
branch when:

- the phenotype evidence is explicit
- the identity limitation is documented clearly
- the row is not misrepresented as strain-resolved

## 5. Habitat-only evidence is not acceptable

### Active project rule

Do not accept an algae row only because it is described as:

- marine
- brackish
- saline environment
- hot spring
- acid spring

These are candidate-discovery clues only.

This is a **project-operational rule**, but it is a necessary one.

### Thesis-ready interpretation

The branch is based on explicit phenotype evidence, not habitat inference.

## 6. Algae salinity labels

### Active project rule

Use:

- `saline-adapted alga` for explicit optimum `0.2 <= NaCl equivalent < 0.85 M`
- `halophilic alga` for explicit optimum `>= 0.85 M`

Keep `moderate` and `extreme` halophile subclasses only when the quantitative
evidence supports them clearly.

These are **project-operational labels**.

### Why this rule is reasonable

The `Asteromonas gracilis` paper shows why explicit quantitative salinity
evidence is feasible and informative in algae:

- [PubMed 16661724](https://pubmed.ncbi.nlm.nih.gov/16661724/)

The `Dunaliella salina SA 134` paper shows why a strain-specific salinity
screen is especially useful when available:

- [PubMed 28808229](https://pubmed.ncbi.nlm.nih.gov/28808229/)

### Thesis-ready interpretation

For algae, the important distinction is not just "salt-tolerant" versus
"salt-associated." The branch should prioritize organisms with an explicit
measured salinity optimum or an experimentally supported growth maximum across
a salinity gradient.

## 7. Algae temperature labels

### Active project rule

Use:

- `high-temperature alga` for explicit optimum growth `>= 40 C`
- `thermophilic alga` for explicit optimum growth `>= 45 C`
- `thermoacidophilic alga` for high-temperature growth plus explicit acidic
  growth context

These are **project-operational labels**.

### Why this rule is reasonable

The verified source pages clearly support the existence of eukaryotic
high-temperature and thermoacidophilic algae:

- [PubMed 34202768](https://pubmed.ncbi.nlm.nih.gov/34202768/)
- [PMC10154280](https://pmc.ncbi.nlm.nih.gov/articles/PMC10154280/)
- [PMC5411153](https://pmc.ncbi.nlm.nih.gov/articles/PMC5411153/)

### Thesis-ready interpretation

The source-backed part is that these organisms occupy high-temperature,
low-pH environments and are biologically extremophilic.

The exact `40 C` and `45 C` branch cutoffs should be described as operational
choices for a reproducible curation workflow.

## 8. Evidence requirement

### Active project rule

Accept a row only when:

- the identity is stable enough to explain
- the phenotype value is explicit enough to defend
- the source is directly citable
- the threshold is satisfied

This is the same curation philosophy used throughout the thesis, adapted for
algae.

## 9. Inclusion decisions

### `include_in_core = yes`

Use `yes` within the algae branch only when:

- identity is defensible
- phenotype evidence is explicit
- the operational threshold is satisfied
- the source is directly citable

### `include_in_core = review`

Use `review` when:

- the candidate is promising
- but identity support is incomplete
- or the full text still needs careful checking
- or the phenotype interpretation is still provisional

### `include_in_core = no`

Use `no` when:

- the threshold is not met
- the evidence is only habitat-based
- the label cannot be defended from the source

## 10. Sequence-linking rule

### Active project rule

Do not start retrieval until the metadata table is stable.

When retrieval starts, use this preference order:

1. exact strain or isolate proteome/genome proteins
2. culture-collection matched genome
3. species-level proteome only if that fallback is documented explicitly

This is a **project-operational retrieval rule**.

## 11. Bottom line for the thesis

The cleanest thesis framing is:

- algae require a separate branch because identity and evidence structure are
  different from the prokaryotic core
- the branch should rely on explicit phenotype evidence, not habitat labels
- source pages justify the biological relevance of saline-adapted, halophilic,
  high-temperature, thermoacidophilic, psychrophilic, and acidophilic algae
- the exact branch thresholds should be declared openly as operational curation
  definitions rather than universal algae standards

## 12. Algae psychrophily label

### Active project rule

Use `psychrophilic_alga` for an alga with explicit optimal growth temperature
`<= 15 C`.

This label is **separate from thermoacidophilic_alga**.

It should not be applied on the basis of habitat description alone (polar,
Arctic, Antarctic, sea-ice). The optimum growth temperature must be explicitly
reported or clearly derivable from a growth-curve experiment.

This is a **project-operational label**.

### Why this rule is reasonable

The `15 C` ceiling mirrors the classical definition of a psychrophile used in
microbiology and environmental genomics: an organism whose optimal growth
temperature is at or below `15 C` and that cannot grow at temperatures above
`20 C` (Morita 1975; Morgan-Kiss et al. 2006).

Two source-backed examples currently in the branch:

- `Fragilariopsis cylindrus` CCMP1102 is a polar marine diatom that grows
  optimally near `4 C`. Its genome paper explicitly reports the psychrophilic
  character of the strain and its adaptation to Antarctic sea ice.
  Source: [Mock et al. 2017, Science 356:1064–1069](https://www.science.org/doi/10.1126/science.aaf4326)

- `Coccomyxa subellipsoidea` C-169 is a green alga isolated from polar
  and high-altitude habitats with documented optimal growth near `5 C`.
  Source: [Blanc et al. 2012, Plant Journal 74:937–946](https://onlinelibrary.wiley.com/doi/10.1111/j.1365-313X.2012.05028.x)

### Thesis-ready interpretation

The source-backed part is that polar and cryophilic algae occupy genuinely
cold-adapted ecological niches and carry molecular signatures of cold
adaptation (cold-active enzymes, membrane fluidity maintenance, antifreeze
proteins).

The exact `<= 15 C` branch cutoff should be described as an operational
choice derived from the accepted microbiology definition.

## 13. Algae acidophily label (non-thermophilic)

### Active project rule

Use `acidophilic_alga` for an alga with explicit optimal growth pH `<= 3.0`
**and without a thermophilic component**.

If an alga satisfies both acid adaptation (`pH <= 3.0`) and high-temperature
growth (`OGT >= 40 C`), it should be classified as `thermoacidophilic_alga`
instead. Do not set both `acidophilic_alga = yes` and
`thermoacidophilic_alga = yes` on the same row.

This label does not apply on habitat inference alone (acidic lake, volcanic
soil). The pH optimum must come from an explicit growth experiment.

This is a **project-operational label**.

### Why this rule is reasonable

The `pH <= 3.0` ceiling is broadly consistent with the extremophile literature
for acidophiles (Baker & Banfield 2003; Seckbach 2007). It excludes organisms
that merely tolerate mildly acidic conditions and restricts the label to
organisms where acid adaptation is a primary ecological trait.

One source-backed example currently in the branch:

- `Chlamydomonas eustigma` NIES-2499 is a strongly acidophilic green alga
  with a documented growth optimum near `pH 1.0`, isolated from highly acidic
  environments. It is mesophilic and therefore classified as acidophilic but
  not thermoacidophilic. Its genome paper documents the genomic basis of acid
  tolerance in a photosynthetic eukaryote.
  Source: [Hirooka et al. 2017, PNAS (PMID 28893987)](https://pubmed.ncbi.nlm.nih.gov/28893987/)

### Thesis-ready interpretation

The source-backed part is that extreme acidophily can occur in photosynthetic
eukaryotes independently of thermophily, and that the molecular strategies for
acid tolerance are distinct from those used by thermoacidophilic Cyanidiales.

The exact `pH <= 3.0` branch cutoff should be described as an operational
choice derived from the standard acidophile definition used across the
extremophile literature.
