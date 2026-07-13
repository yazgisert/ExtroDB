# Step 3: Screen TEMPURA Thermophile Candidates

You already have the filtered TEMPURA result with `Topt_average >= 45 C`.

At this point, do **not** go to sequence collection yet.

The job now is to convert the TEMPURA export into a **screening list** and then
promote the strongest entries into the core taxon table.

## Goal of this step

Produce:

- `10 to 20` TEMPURA candidates in the screening log
- `~10` strong thermophile rows ready for BacDive validation

## Use these TEMPURA columns first

From the TEMPURA export, focus on:

- `Genus_and_species`
- `Taxonomy_ID`
- `Strain`
- `Topt_average(°C)` or equivalent column name
- optionally `Topt_low(°C)` and `Topt_high(°C)`

Ignore the rest for the first pass unless needed later.

## Selection rule for the first shortlist

Prefer rows that have all of the following:

- non-empty `Genus_and_species`
- non-empty `Strain`
- non-empty `Taxonomy_ID`
- numeric `Topt_average`

If you have too many rows, prioritize:

1. rows with clear strain names
2. rows with high `Topt_average`
3. rows with species-level names that look stable

## How to classify TEMPURA rows

Based on `Topt_average`:

- `45 <= Topt_average < 80` -> `thermophile`
- `Topt_average >= 80` -> `hyperthermophile`

For this TEMPURA stage, use:

- `evidence_level_candidate = curated_database`
- `follow_up_action = validate_in_BacDive`

## What to copy into the screening log

For each selected TEMPURA row, fill:

- `source_type = TEMPURA`
- `citation_or_resource_url = https://togodb.org/db/tempura`
- `focus_category = thermophile`
- `taxon_name_candidate = Genus_and_species`
- `strain_candidate = Strain`
- `key_evidence = TEMPURA Topt_average = X C`
- `extracted_ogt_celsius = Topt_average`
- `extracted_nacl_optimum_M =` leave empty
- `extracted_label = thermophile` or `hyperthermophile`
- `evidence_level_candidate = curated_database`
- `follow_up_action = validate_in_BacDive`
- `notes =` optional comments

## What not to do yet

Do not do these yet:

- decide final inclusion solely from TEMPURA
- collect protein sequences
- add AFESM candidates
- clean every taxonomy issue in full detail

## When a row is good enough to validate in BacDive

A TEMPURA row is ready for BacDive if:

- the species name is readable
- the strain name is not empty
- the temperature value is numeric

That is enough for the next step.

## First-pass target

Your first-pass target is:

- pick `20` TEMPURA rows
- then BacDive-check the best `10`

This keeps the workload realistic and lets you test the pipeline early.
