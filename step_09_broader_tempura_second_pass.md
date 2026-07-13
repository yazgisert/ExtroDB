# Step 9: Broader TEMPURA Second Pass

This step does not weaken the thermophile inclusion criteria.

The thermophile rules stay the same:

- `thermophile`: `OGT >= 45 C`
- `hyperthermophile`: `OGT >= 80 C`
- only strong strain-aware evidence becomes `include_in_core = yes`

What changes in this step is only the *candidate discovery width*.

## Why a broader second pass is needed

The first local TEMPURA plus BacDive pass has already been largely exhausted:

- `40` TEMPURA candidates were screened
- `33` became the main core
- `3` more were promoted into a small round-2 batch
- the remaining unresolved rows are edge cases

So if you want many more thermophile candidates, the next move is:

- go back to a broader TEMPURA export
- generate a larger second-pass screening list
- validate those rows in BacDive under the same strict rules

## What input is needed

You need a broader TEMPURA export file that still includes at least:

- `Genus_and_species`
- `Strain`
- `Taxonomy_ID`
- `Topt_average`

The export can be TSV or CSV.

## Script

Use:

- `prepare_tempura_second_pass.py`

This script:

- reads a broader TEMPURA export
- keeps only rows with usable species, strain, and `Topt_average`
- keeps only rows with `Topt_average >= 45`
- excludes rows already present in `paper_screening_log_tempura_bacdive.tsv`
- labels each remaining row as `thermophile` or `hyperthermophile`
- writes a new screening log for the next BacDive validation pass

## Example command

```bash
python3 /Users/yazgisert/Documents/Playground/thesis_data_collection/prepare_tempura_second_pass.py \
  /path/to/tempura_export.tsv \
  --output /Users/yazgisert/Documents/Playground/thesis_data_collection/paper_screening_log_tempura_broad_second_pass.tsv \
  --limit 50
```

## Recommended workflow

1. Export a broader TEMPURA candidate table.
2. Run `prepare_tempura_second_pass.py`.
3. Review the resulting broad second-pass screening log.
4. Validate the best rows in BacDive.
5. Promote only strong rows into a thermophile expansion core table.
6. Then do GTDB-first sequence linking and protein download.

## Important interpretation

This is a broader *screening* step, not a broader *acceptance* step.

That distinction is the key methodological point:

- broader candidate pool
- same strict inclusion logic
