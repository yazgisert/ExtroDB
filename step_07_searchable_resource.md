# Step 7: Searchable Dataset Layer

The next step is not to replace the curated tables.
It is to add a thin searchable layer on top of the current prokaryotic core.

That searchable layer should preserve the logic already used in the thesis:

- keep the prokaryotic thermophile and halophile core as the active dataset
- keep phenotype evidence separate from sequence retrieval status
- keep raw FASTA files, manifests, and merged FASTAs linked back to the
  curation table

## Why this is useful

This solves two needs at once:

- it gives one place to filter the dataset by phenotype, class, evidence level,
  or inclusion decision
- it creates the first minimal version of a resource that could later be shared
  with other researchers

## What gets combined

The searchable index joins the current active sources:

- thermophile core curation table
- thermophile sequence retrieval tracker
- thermophile manifest
- halophile core curation table
- halophile sequence retrieval tracker
- halophile manifest

Each output row still represents one curated organism or strain-level entry.
This is still an organism-level proteome resource, not a protein-level label
resource.

## Output idea

The searchable layer should answer questions like:

- show all downloaded extreme halophiles
- show all thermophiles with OGT at or above a threshold
- show rows with missing proteomes

So the output needs to keep:

- biological labels
- evidence fields
- inclusion decision
- sequence retrieval status
- FASTA and manifest pointers

## Current implementation

The current builder is:

- `build_searchable_extremophile_index.py`

It writes:

- `searchable_resource/extremophile_master_index.tsv`
- `searchable_resource/extremophile_master_summary.tsv`

The query helper is:

- `query_extremophile_master_index.py`

## Recommended workflow

1. Update the curation tables or retrieval trackers first.
2. Rebuild manifests if proteome downloads changed.
3. Rebuild the searchable index.
4. Query the index instead of manually editing the source tables.

The current searchable index intentionally excludes the eukaryotic algae
branch. Algae will be restarted later under separate criteria.

## Example commands

Build the index:

```bash
python3 /Users/yazgisert/Documents/Playground/thesis_data_collection/build_searchable_extremophile_index.py
```

Count hyperthermophile rows:

```bash
python3 /Users/yazgisert/Documents/Playground/thesis_data_collection/query_extremophile_master_index.py \
  --where dataset_branch=thermophile_core \
  --where hyperthermophile=yes \
  --count
```

Show downloaded extreme halophiles:

```bash
python3 /Users/yazgisert/Documents/Playground/thesis_data_collection/query_extremophile_master_index.py \
  --where dataset_branch=halophile_core \
  --where halophile_class=extreme \
  --where proteome_status=downloaded \
  --columns taxon_name,strain,protein_count,local_fasta_path
```

Show accepted rows with missing proteomes:

```bash
python3 /Users/yazgisert/Documents/Playground/thesis_data_collection/query_extremophile_master_index.py \
  --where proteome_status=missing \
  --columns dataset_branch,taxon_name,strain,retrieval_notes
```
