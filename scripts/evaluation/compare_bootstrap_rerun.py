# Checks the rerun of the 5-class and 12-class models with 1000 bootstrap iterations (instead of 30): the point values
# (test metrics and per-class metrics) must be identical to the existing results, so only the confidence intervals change.
# Prints the old and new macro-F1 CIs. With --apply, the rerun files replace the existing ones, but only if all match.
# Run from the repository root: python scripts/evaluation/compare_bootstrap_rerun.py [--apply]

import argparse
import json
import shutil
from pathlib import Path

OLD_DIR = Path("results/training_results")
NEW_DIR = Path("results/training_results_rerun")
SPLITS  = ["split1_30pct", "split1_20pct", "split2_30pct", "split2_20pct"]
SUFFIX  = {"ProtT5": "", "ESM2 650M": "_esm2", "ESM2 8M": "_esm2_8m", "Composition": "_onehot"}


def same(a, b, tol=1e-9):
    if isinstance(a, dict):
        return a.keys() == b.keys() and all(same(a[k], b[k], tol) for k in a)
    return abs(a - b) <= tol


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--apply", action="store_true")
    apply = parser.parse_args().apply

    all_match, files = True, []
    print(f"{'split':<14} {'embedder':<12} {'macro-F1':>9}  {'CI (30 iterations)':>18}  {'CI (1000 iterations)':>20}  points")
    for split in SPLITS:
        for embedder, suffix in SUFFIX.items():
            name = f"{split}{suffix}_result.json"
            if not (NEW_DIR / name).exists():
                print(f"{split:<14} {embedder:<12} missing in {NEW_DIR}")
                all_match = False
                continue
            old, new = json.load(open(OLD_DIR / name)), json.load(open(NEW_DIR / name))
            match = same(old["test_metrics"], new["test_metrics"]) and same(old["per_class_metrics"], new["per_class_metrics"])
            all_match &= match
            files.append(name)
            o, n = old["test_metrics_bootstrap"]["macro-f1_score"], new["test_metrics_bootstrap"]["macro-f1_score"]
            print(f"{split:<14} {embedder:<12} {new['test_metrics']['macro-f1_score']:>9.3f}  "
                  f"{o['lower']:>8.3f}–{o['upper']:<9.3f}  {n['lower']:>9.3f}–{n['upper']:<10.3f}  "
                  f"{'identical' if match else 'DIFFERENT'}")

    if not all_match:
        print("\nNot all point values match; the existing results are kept.")
    elif apply:
        for name in files:
            shutil.copy(NEW_DIR / name, OLD_DIR / name)
        print(f"\nAll point values identical; {len(files)} files copied to {OLD_DIR}.")
    else:
        print("\nAll point values identical. Run again with --apply to replace the existing results.")
