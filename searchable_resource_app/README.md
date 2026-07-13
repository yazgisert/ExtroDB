# Searchable Resource App

This local app adds a minimal browser interface on top of the curated
`extremophile_master_index.tsv` resource.

The master index is proteome-deduplicated so each downloaded proteome is counted
once even when the same strain belongs to multiple extremophile labels. Branch-
level provenance is preserved separately in
`searchable_resource/extremophile_membership_index.tsv`.

## Run

```bash
python3 /Users/yazgisert/Documents/Playground/thesis_data_collection/searchable_resource_webapp.py
```

Then open:

```text
http://127.0.0.1:8123
```

## What it does

- loads the current deduplicated searchable resource directly from the local TSV
- filters rows by branch, phenotype, proteome status, thermal class, and
  numeric thresholds
- searches across taxon names, strain labels, accession IDs, and notes
- shows a detail panel with evidence, retrieval notes, and local file links
- exports the current filtered selection to a TSV file

## Scope

This is intentionally a small local research tool.
It reflects the current organism-level, proteome-deduplicated index, not a
protein-level annotation database.
