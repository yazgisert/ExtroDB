# Sources for Thresholds and Curation Criteria

This note separates:

1. literature-backed biological definitions
2. project-specific operational rules

That distinction is important in the thesis because not every rule in the dataset is a universal community standard. Some are standard microbiology definitions, while others are explicit curation decisions chosen to make the dataset rigorous and reproducible.

## 0. Current project state

At the current stage of the project, the thesis workflow has three clearly
separated layers:

1. an active prokaryotic core
2. a searchable resource built from that prokaryotic core
3. a reset algae branch that will restart later under separate rules

### Active prokaryotic core

The current active thesis dataset includes:

- prokaryotic psychrophiles
- prokaryotic thermophiles
- prokaryotic halophiles
- prokaryotic acidophiles
- prokaryotic moderate acidophiles
- prokaryotic acidotolerants

These are the phenotype branches currently supported by the working curation
and sequence pipeline.

### Searchable resource layer

The searchable resource is a convenience layer built on top of the active
prokaryotic core.
It joins:

- curation tables
- sequence retrieval trackers
- manifest tables
- merged FASTA outputs

It does not change the biological inclusion logic.
It only makes the active prokaryotic core easier to inspect and query.

### Algae branch status

The earlier exploratory algae rows have now been archived and the algae TSVs
have been reset to empty starting tables.

That means:

- algae are not part of the active searchable core
- earlier algae examples should be treated as historical exploratory material
- the restarted algae branch should begin again from explicit discovery and
  strain-reconciliation steps

### Algae restart source stack

The planned restart order for algae is now:

1. `AlgaeBase` for candidate discovery, accepted names, and synonym cleanup
2. `UTEX`, `CCAP`, and `NCMA` for strain-level identity reconciliation
3. primary literature for explicit physiological evidence
4. review literature only as supporting context

This is important for the thesis because it gives a cleaner justification than
the older exploratory algae branch, which was more literature-first and less
systematically strain-anchored.

## 1. Prokaryotic core: thermophiles

### Active rule in this project

- `thermophile`: `OGT >= 45 C`
- `hyperthermophile`: `OGT >= 80 C`

### Literature support

1. Frontiers review:
   [Extremely thermophilic microorganisms as metabolic engineering platforms for production of fuels and industrial chemicals](https://www.frontiersin.org/journals/microbiology/articles/10.3389/fmicb.2015.01209/full)

   This review states that organisms with optimum temperature above `45 C` are thermophiles, and that hyperthermophiles grow optimally at `80 C` and above.

2. Review on hyperthermophilic enzymes:
   [Hyperthermophilic Enzymes: Sources, Uses, and Molecular Mechanisms for Thermostability](https://pmc.ncbi.nlm.nih.gov/articles/PMC99017/)

   This review defines hyperthermophiles as bacteria and archaea with optimal growth temperatures above `80 C`.

3. ThermoBase paper:
   [ThermoBase: A database of the phylogeny and physiology of thermophilic and hyperthermophilic organisms](https://pmc.ncbi.nlm.nih.gov/articles/PMC9089862/)

   This paper is useful methodologically because it explains that there is no fully universal cutoff in the literature, with thermophile definitions appearing at `40 C`, `45 C`, or `50 C`, and hyperthermophile definitions at `70 C`, `75 C`, or `80 C`.

### Thesis-ready interpretation

You can defend the `45 C` threshold as a conservative, widely used cutoff for thermophily, and the `80 C` threshold as a standard cutoff for hyperthermophily. The ThermoBase paper is especially useful because it shows that threshold ambiguity exists in the literature, which justifies stating your chosen cutoff explicitly as an operational definition.

## 1A. Psychrophile branch

### Branch rule in this project

- `psychrophile`: `OGT <= 15 C`
- exclude `psychrotolerant` or `psychrotrophic` organisms from the strict
  psychrophile branch when the optimum is above `15 C`, even if low-temperature
  growth is reported

### Literature support

1. Classical review:
   [Psychrophilic bacteria](https://doi.org/10.1128/MMBR.39.2.144-167.1975)

   Morita is one of the standard sources for the distinction between true
   psychrophiles and organisms that merely tolerate low temperature.

2. Review on cold adaptation:
   [Psychrophilic microorganisms: challenges for life](https://doi.org/10.1038/sj.embor.7400662)

   This review is useful for thesis framing because it summarizes the commonly
   used cold-side definition and emphasizes the distinction between strict
   psychrophiles and psychrotolerant or psychrotrophic organisms.

3. Strain-level example of a defensible psychrophile:
   [Alpinimonas psychrophila gen. nov., sp. nov., an actinobacterium of the family Microbacteriaceae isolated from alpine glacier cryoconite](https://doi.org/10.1099/ijs.0.036160-0)

   This species description reports good growth over `1-15 C`, which is a
   clean example of the kind of cold-adapted organism that fits a strict
   psychrophile branch.

### Thesis-ready interpretation

The cleanest project rule is `OGT <= 15 C` for a strict psychrophile label.
That aligns with the widely used classical cold-side definition and also fits
the kind of strain-level examples that can be defended in curation.

In the current psychrophile workflow, the practical source stack now includes:

- `TEMPURA` for exact-strain OGT rows when the local snapshot links a genome
  assembly directly
- `BacDive` for exact type-strain culture-temperature entries used in the
  expansion round
- primary literature for explicit growth ranges or optimum temperatures when a
  matched exact-strain assembly is available

Because the current thesis tables store optimum growth temperature much more
consistently than maximum growth temperature, the most practical operational
choice is to use `OGT <= 15 C` as the main threshold. Organisms that grow at
low temperature but have optima above `15 C` should be treated as
`psychrotolerant` or `psychrotrophic`, not merged into the strict
psychrophile branch.

For branch design, psychrophiles fit best as a separate cold-temperature
branch, or as a sibling branch under a broader `temperature-related`
umbrella. They should not be folded into a branch that remains explicitly
named `thermophile`, because that would mix opposite temperature-adaptation
regimes under a misleading label.

The branch is now active in the project as a separate workflow with a first
curated pass and a prepared sequence-retrieval tracker.

## 2. Prokaryotic core: halophiles

### Active rule in this project

- `halophile`: `NaCl optimum > 0.2 M`
- `slight halophile`: `0.2 <= NaCl optimum < 0.85 M`
- `moderate halophile`: `0.85 <= NaCl optimum < 3.4 M`
- `extreme halophile`: `NaCl optimum >= 3.4 M`

### Literature support

1. Review with the exact class boundaries used in this project:
   [The Methods of Digging for “Gold” within the Salt: Characterization of Halophilic Prokaryotes and Identification of Their Valuable Biological Products Using Sequencing and Genome Mining Tools](https://www.mdpi.com/2073-4425/12/11/1756)

   This review gives:
   - slight halophiles: `0.2-0.85 M`
   - moderate halophiles: `0.85-3.4 M`
   - extreme halophiles: `3.4-5.1 M`

   It also distinguishes halophiles from halotolerant organisms by salt requirement for optimal growth.

2. Background review on high-salt life:
   [Microbial life at high salt concentrations: phylogenetic and metabolic diversity](https://aquaticbiosystems.biomedcentral.com/articles/10.1186/1746-1448-4-2)

   This paper is useful for thesis context because it reviews classic halophile terminology and shows that halophile classification has a long literature history rather than being an ad hoc choice.

### Thesis-ready interpretation

You can defend `NaCl optimum > 0.2 M` as the operational boundary that excludes halotolerant organisms and retains organisms whose optimal growth is genuinely saline-adapted. The subclass boundaries are directly aligned with widely used halophile class ranges in the literature.

## 3. Prokaryotic core: acid-related pH labels

### Active rule in this project

- `acidophile`: optimum growth `pH <= 3.0`
- `moderate acidophile`: optimum growth `3.0 < pH <= 5.0`
- `acidotolerant`: optimum growth `5.0 < pH <= 7.5`, but able to grow at lower pH
- optional project-specific subclass:
  - `hyperacidophile`: optimum growth `pH <= 1.0`

### Literature support

1. Foundational review:
   [Life in acid: pH homeostasis in acidophiles](https://pubmed.ncbi.nlm.nih.gov/17331729/)

   This review gives a clean thesis-usable definition by stating that microorganisms with optimum growth below `pH 3` are termed acidophiles.

2. Editorial overview:
   [Editorial: Recent Advances in Acidophile Microbiology: Fundamentals and Applications](https://www.frontiersin.org/article/10.3389/fmicb.2017.00428/full)

   This source gives the full three-part distinction that is most useful for curation:
   - extreme acidophiles: optimum growth at or below `pH 3.0`
   - moderate acidophiles: optimum growth between `pH 3.0` and `5.0`
   - acid-tolerant organisms: optimum growth above `pH 5.0`, but able to grow reasonably well at lower pH

3. FEMS review:
   [Geomicrobiology of extremely acidic subsurface environments](https://pubmed.ncbi.nlm.nih.gov/22224750/)

   This review supports the same framework and is especially useful because it distinguishes clearly between truly acidophilic organisms and acid-tolerant organisms that remain active in acidic environments without having an acidic pH optimum.

4. Example of a very low-pH specialist:
   [Picrophilus gen. nov., fam. nov.: a novel aerobic, heterotrophic, thermoacidophilic genus and family comprising archaea capable of growth around pH 0](https://pmc.ncbi.nlm.nih.gov/articles/PMC177581/)

   This primary paper is useful if you want an explicit example for an optional very-low-pH subclass, because it reports a `pH` optimum of `0.7` for `Picrophilus`.

5. Primary example for the acidotolerant branch:
   [Bryobacter aggregatus gen. nov., sp. nov., a peat-inhabiting, aerobic chemo-organotroph from subdivision 3 of the Acidobacteria](https://pubmed.ncbi.nlm.nih.gov/19651730/)

   This paper is especially useful for the revised acidotolerant branch because it explicitly describes the strains as `acidotolerant` and reports growth at `pH 4.5-7.2` with an optimum at `pH 5.5-6.5`. That is a clean primary example of the exact phenotype category you now want to include.

### Thesis-ready interpretation

These sources support using `pH <= 3.0` as the main operational boundary for acidophily if acidophiles are added to the prokaryotic core. They also show that it would be too strict to classify every organism with optimum growth above `pH 3.0` as merely acidotolerant, because many sources recognize `moderate acidophiles` with optima between `pH 3` and `5`.

For the revised acidotolerant branch, the literature-backed part is the lower boundary:
- reviews distinguish acid-tolerant organisms from acidophiles by placing acid-tolerant organisms above `pH 5.0`
- primary papers such as `Bryobacter aggregatus` show concrete acidotolerant taxa with optima in the `pH 5.5-6.5` range

For this reason, the cleanest thesis framing is:

- use `acidophile` as the core label for organisms with optimum growth at `pH <= 3.0`
- optionally track `moderate acidophile` separately if you want broader coverage later
- reserve `acidotolerant` for organisms whose optimum growth is above `pH 5.0`
- in this project, cap the acidotolerant branch at `pH <= 7.5` as an explicit operational choice rather than a canonical literature cutoff, so it remains focused on organisms with still-acid-shifted growth optima rather than broad near-neutral generalists
- treat `pH <= 1.0` as an optional project-specific subclass such as `hyperacidophile`, rather than calling it the only `extreme acidophile` class, because many reviews already use `extreme acidophile` for the broader `pH <= 3.0` group

## 4. Prokaryotic core: alkaliphiles

### Active rule in this project

- `alkaliphile`: optimum growth `pH >= 9.0`
- optional project-specific subclass:
  - `extreme alkaliphile`: optimum growth `pH >= 10.0`
- `alkalitolerant`: optimum growth `7.5 < pH < 9.0`, but with documented
  growth at `pH >= 9.0`

### Literature support

1. Review with a directly usable definition:
   [Alkaliphiles](https://www.jstage.jst.go.jp/article/pjab/80/4/80_4_166/_article)

   In the abstract, Horikoshi defines alkaliphiles as microorganisms that
   grow optimally or very well at `pH > 9` and that cannot grow, or grow only
   slowly, near `pH 6.5`.

2. Widely cited review summarized in later secondary sources:
   `Horikoshi K. Alkaliphiles: Some applications of their products for
   biotechnology. Microbiology and Molecular Biology Reviews. 1999.`

   This review is consistently cited as describing alkaliphiles as organisms
   adapted to alkaline environments with growth optima around `pH 10`, which
   is useful context for a stricter subclass.

3. Primary obligate-alkaliphile example:
   `Yumoto I, Yamaga S, Sogabe Y, Nodasaka Y, Matsuyama H. Bacillus
   krulwichiae sp. nov., a halotolerant obligate alkaliphile that utilizes
   benzoate and m-hydroxybenzoate. International Journal of Systematic and
   Evolutionary Microbiology. 2003.`

   This is useful as a strain-level example of a taxon explicitly described as
   an `obligate alkaliphile`, not just alkali-tolerant.

4. Second primary obligate-alkaliphile example:
   `Yumoto I, Hirota K, Kawahara T, Nodasaka Y, Okuyama H, Matsuyama H,
   Yokota Y, Hoshino T, Nakajima K. Oceanobacillus oncorhynchi sp. nov., a
   halotolerant obligate alkaliphile isolated from the skin of a rainbow trout
   (Oncorhynchus mykiss), and emended description of the genus Oceanobacillus.
   International Journal of Systematic and Evolutionary Microbiology. 2005.`

   This gives another thesis-usable example that the obligate alkaliphile
   concept is used explicitly at the strain-species description level.

5. Primary alkali-tolerant example:
   `Abildgaard L, Nielsen MB, Kjeldsen KU, Ingvorsen K. Desulfovibrio
   alkalitolerans sp. nov., a novel alkalitolerant, sulphate-reducing
   bacterium isolated from district heating water. International Journal of
   Systematic and Evolutionary Microbiology. 2006.`

   This is useful because it shows that `alkalitolerant` is also used as a
   distinct phenotype label rather than being interchangeable with
   `alkaliphile`.

### Thesis-ready interpretation

The most defensible literature-backed boundary is the main `alkaliphile`
threshold:

- use `pH >= 9.0` as the hard entry criterion for alkaliphily

That part is directly supported by Horikoshi's review definition.

The following parts should be presented as explicit project-level operational
choices inferred from the broader literature rather than as universal canonical
cutoffs:

- `extreme alkaliphile = pH >= 10.0`
- `alkalitolerant = 7.5 < pH < 9.0` with documented growth at `pH >= 9.0`

That framing is strong because:

- many classic alkaliphile examples cluster around `pH 10`
- the literature clearly distinguishes true alkaliphiles from merely
  alkali-tolerant organisms
- but it does not give one universally fixed numerical cutoff for every
  possible alkalitolerant subclass

So the clean thesis wording is:

- `alkaliphile` is literature-backed at `pH >= 9.0`
- `extreme alkaliphile` is a useful project-specific subclass for stronger
  alkaline specialists
- `alkalitolerant` should be reserved for organisms with sub-alkaliphilic
  optima that still show documented alkaline growth
- organisms whose optimum remains neutral or only weakly alkaline should stay
  outside the alkaliphile core

## 5. Prokaryotic core: evidence and inclusion rules

### Active rule in this project

- only rows with explicit phenotype evidence and defensible taxon mapping become `include_in_core = yes`
- ambiguous strain mapping, borderline values, or species-level-only cases stay `review`
- GTDB is used first for sequence linking, with NCBI Assembly fallback when needed

### Literature support

1. ThermoBase methodological point:
   [ThermoBase](https://pmc.ncbi.nlm.nih.gov/articles/PMC9089862/)

   This paper explicitly notes that thermophile and hyperthermophile labels are hard to assign reliably without going back to the original literature on optimal growth temperature.

2. Dunaliella ecology review:
   [The ecology of Dunaliella in high-salt environments](https://jbiolres.biomedcentral.com/counter/pdf/10.1186/s40709-014-0023-y.pdf)

   This review is very useful for your general argument that habitat alone is not enough to infer phenotype. It shows that Dunaliella species often live in environments with salt concentrations much higher than their laboratory growth optima. In other words, occurrence in a hypersaline habitat does not automatically tell you the organism's optimum salinity.

### Thesis-ready interpretation

These papers support the logic of demanding explicit physiological evidence rather than relying on habitat words such as `hot spring`, `saltern`, or `salt lake`. That exact requirement is still your project rule, but it is literature-informed and methodologically defensible.

## 6. Exploratory algae branch: what is literature-backed vs operational

The algae branch is different from the prokaryotic core.

For algae, the framework in `Operational Definitions 2` is best described as:

- literature-informed
- exploratory
- operational rather than canonical

That means some parts come from established biology, while others are your explicit curation choices for this thesis.

## 7. Algae: salinity labels

### Active exploratory rule

- `saline-adapted alga`: explicit salinity optimum `0.2 <= x < 0.85 M`
- `halophilic alga`: explicit salinity optimum `>= 0.85 M`
- `AlgaeBase` is used for discovery and synonym control, not as a phenotype
  source by itself
- `UTEX`, `CCAP`, and `NCMA` are used as identity-supporting collection sources
  when strain-level records exist

### Literature support

1. Same halophile-class source as above:
   [Methods of Digging for “Gold” within the Salt](https://www.mdpi.com/2073-4425/12/11/1756)

   The boundaries themselves are standard halophile class boundaries and are not restricted to prokaryotes.

2. Dunaliella review:
   [The ecology of Dunaliella in high-salt environments](https://jbiolres.biomedcentral.com/counter/pdf/10.1186/s40709-014-0023-y.pdf)

   This review reports that:
   - `D. viridis` grows optimally at about `60-90 g/L NaCl`
   - `D. salina` grows optimally at about `100-150 g/L NaCl`

   It also emphasizes that these optima are often much lower than the salinity of the habitats where these algae are found.

3. Dunaliella maritima paper:
   [Cloning and Characterization of Two Putative P-Type ATPases from the Marine Microalga Dunaliella maritima](https://www.mdpi.com/2223-7747/10/12/2667)

   This paper summarizes that some Dunaliella species have growth optima around `2-4% NaCl`, while hypersaline or halophilic Dunaliella species have optima above `6% NaCl`.

4. Experimental paper on `D. salina`:
   [Optimization of Salt Concentrations for a Higher Carotenoid Production in Dunaliella salina](https://pubmed.ncbi.nlm.nih.gov/27020189/)

   This paper reports maximal growth for one `D. salina` study strain at `1.5-3.0 M NaCl`.

### Thesis-ready interpretation

The algae salinity labels are defensible as a transfer of the standard halophile salinity framework into a separate eukaryotic exploratory dataset. The important caveat is that this is an operational classification scheme for your algae branch, not a universal formal standard for algal ecology.

The stronger methodological point now is that algae discovery and strain
reconciliation should be separated.
`AlgaeBase` helps define the candidate name space, while `UTEX`, `CCAP`, and
`NCMA` help determine whether a candidate can be defended at strain or isolate
level before sequence retrieval begins.

## 7. Algae: temperature labels

### Active exploratory rule

- `high-temperature alga`: explicit growth optimum `>= 40 C`
- `thermophilic alga`: explicit growth optimum `>= 45 C`
- `thermoacidophilic alga`: temperature rule plus explicit acidic growth context

### Literature support

1. Classic physiological study:
   [The Physiological Ecology of Cyanidium caldarium](https://www.microbiologyresearch.org/content/journal/micro/10.1099/00221287-67-1-17)

   This study reports that `Cyanidium caldarium` was found in habitats with temperatures above `40 C`, had an upper limit around `55-56 C`, and had an optimum growth temperature of `45 C`.

2. Applied and Environmental Microbiology paper:
   [Biogeographic and Phylogenetic Diversity of Thermoacidophilic Cyanidiales in Yellowstone National Park, Japan, and New Zealand](https://journals.asm.org/doi/abs/10.1128/aem.02741-07)

   This paper describes Cyanidiales as phototrophs living in environments that combine low pH with moderately high temperatures of about `40-56 C`.

3. Scientific Reports paper:
   [Efficient open cultivation of cyanidialean red algae in acidified seawater](https://www.nature.com/articles/s41598-020-70398-z)

   This paper describes Cyanidiales as dominating sulfuric acidic hot springs worldwide at about `pH 0.05-5.0` and `35-56 C`.

4. BMC Ecology and Evolution paper:
   [Establishment of endolithic populations of extremophilic Cyanidiales](https://bmcecolevol.biomedcentral.com/articles/10.1186/1471-2148-6-78)

   This paper describes Cyanidiales as eukaryotic extremophiles growing at high temperature around `50-55 C` and very low pH around `0.5-3.0`.

### Thesis-ready interpretation

The `>= 45 C` algae thermophily rule is easy to defend because it aligns with the prokaryotic thermophile cutoff and is also consistent with the classic `Cyanidium caldarium` optimum of `45 C`.

The `>= 40 C` `high-temperature alga` label is more of a deliberate exploratory category. It is defensible because thermoacidophilic red algae are repeatedly described in the literature as inhabiting about `35-56 C` or `40-56 C` systems, but it should be presented as your own broader screening label rather than a universal textbook definition.

## 8. Algae: why species-level evidence may be acceptable

### Active exploratory rule

- culture-collection strains are preferred
- named isolates are acceptable
- species-level entries may be acceptable if phenotype evidence is explicit and limitations are documented

For the restarted algae branch, this should now be interpreted using a concrete
source ladder:

1. collection-backed strain or isolate if available
2. named study strain
3. species-level fallback only when strain resolution is not realistically
   recoverable

### Source status

This is mainly a project-specific methodological choice, not a universally fixed biological threshold from one paper.

### Thesis-ready interpretation

You can defend this by saying that the algae branch is exploratory and that the eukaryotic literature often provides strong phenotype information without the same strain-standardization quality available for prokaryotic BacDive-style curation. Therefore, species-level acceptance was allowed only in the exploratory algae branch, never in the thesis-grade prokaryotic core, and only when the phenotype evidence itself was explicit.

## 9. Algae: why habitat-only language is insufficient

### Active exploratory rule

- habitat terms alone such as `salt lake`, `saltern`, `hot spring`, or `acid spring` are not enough
- explicit growth, optimum, range, or experimentally supported phenotype evidence is preferred

### Literature support

1. Dunaliella ecology review:
   [The ecology of Dunaliella in high-salt environments](https://jbiolres.biomedcentral.com/counter/pdf/10.1186/s40709-014-0023-y.pdf)

   This is the strongest source for this point because it shows directly that environmental salinity and growth optimum are not the same thing.

2. ThermoBase:
   [ThermoBase](https://pmc.ncbi.nlm.nih.gov/articles/PMC9089862/)

   This paper supports the same logic on the temperature side: reliable thermophile assignment depends on actual physiological growth data rather than loose labeling.

### Thesis-ready interpretation

This rule is strongly defensible and should be stated as a curation-quality requirement rather than as a biological claim.

## 10. Suggested thesis wording

You can say something close to this:

`For the prokaryotic core dataset, thermophily was operationally defined as an optimal growth temperature >= 45 C and hyperthermophily as >= 80 C, consistent with widely used definitions in the extremophile literature, while halophily was defined using NaCl optimum-based classes, with > 0.2 M required for inclusion as halophilic. Because published thresholds vary somewhat across sources, all cutoffs were stated explicitly and only rows supported by direct physiological evidence and defensible taxon mapping were accepted into the core dataset.`

For algae:

`Because the eukaryotic algae literature is less consistently strain-resolved than the prokaryotic literature, algae were treated as a separate exploratory dataset with literature-informed but explicitly operational criteria. Species-level acceptance was therefore allowed only in the exploratory algae branch and only when phenotype evidence was explicit, while habitat-only descriptions were not considered sufficient for classification.`

You can now also add a source-structure sentence such as:

`For the restarted algae branch, candidate discovery was separated from phenotype validation: AlgaeBase was used to normalize names and synonyms, culture-collection resources such as UTEX, CCAP, and NCMA were used to improve strain-level identity reconciliation where possible, and phenotype labels were assigned only from explicit physiological evidence rather than from habitat language alone.`
