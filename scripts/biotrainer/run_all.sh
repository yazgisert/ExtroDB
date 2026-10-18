#!/usr/bin/env bash
# Trains every prepared predictor with biotrainer (one run per config.yml); finished runs are skipped.
# Run from the repository root in the biotrainer env: bash scripts/biotrainer/run_all.sh
set -e
for config in biotrainer_runs/*/*/config.yml; do
    run_dir=$(dirname "$config")
    if [ -d "$run_dir/output" ]; then
        echo "skip $run_dir (output exists)"
        continue
    fi
    echo "=== $run_dir"
    biotrainer train --config "$config" > "$run_dir/train.log" 2>&1
    echo "    done, log in $run_dir/train.log"
done
