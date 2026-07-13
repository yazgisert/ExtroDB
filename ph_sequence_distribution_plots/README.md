# pH Sequence Distribution Plots

This folder contains simple sequence-count plots for the current pH-related
datasets.

## What is plotted

- `acid_side_sequence_distribution.svg`
  - built from:
    - `acidophile_manifests/sequences_manifest_combined_acidophiles_36.tsv`
    - `acidotolerant_manifests/sequences_manifest_acidotolerant_10.tsv`
- `alkaline_side_sequence_distribution.svg`
  - built from:
    - `alkaliphile_manifests/sequences_manifest_alkaliphile_16.tsv`

The bars show **protein sequence counts** aggregated by the stored
`ph_optimum` value in the manifest files.

## Files

- `plot_ph_sequence_distribution.py`: plotting script
- `*_sequence_distribution.svg`: plots
- `*_sequence_counts_by_ph.csv`: exact pH totals
- `*_sequence_counts_by_interval.csv`: totals for the pH intervals shown in the
  plot backgrounds

## How to run

Use the bundled Python runtime from the Codex desktop app:

```bash
/Users/yazgisert/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/bin/python3 \
  /Users/yazgisert/Documents/Playground/thesis_data_collection/ph_sequence_distribution_plots/plot_ph_sequence_distribution.py
```

The script uses only the Python standard library and writes SVG files directly,
so it does not need `matplotlib`.
