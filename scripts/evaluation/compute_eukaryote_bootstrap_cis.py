# Adds 95% bootstrap CIs to the binarized 5-class algae results. Recall is a fraction of n proteins, so the
# predictions can be rebuilt from n and recall; no model run is needed.

import json
from pathlib import Path

import numpy as np

N_BOOTSTRAP = 1000
SEED        = 42
CI          = 0.95
RNG         = np.random.default_rng(SEED)

RESULTS_DIR = Path("results/training_results")

BINARIZED_JSON  = RESULTS_DIR / "binarized_fiveclass_eukaryote_eval_result.json"


def bootstrap_recall_ci(n_proteins: int, recall: float) -> tuple[float, float]:
    n_correct = round(n_proteins * recall)
    preds     = np.array([1] * n_correct + [0] * (n_proteins - n_correct))
    recalls   = np.array([
        preds[RNG.integers(0, n_proteins, size=n_proteins)].mean()
        for _ in range(N_BOOTSTRAP)
    ])
    lo = (1 - CI) / 2 * 100
    hi = (1 - (1 - CI) / 2) * 100
    return float(np.percentile(recalls, lo)), float(np.percentile(recalls, hi))


def add_class_cis(by_class: dict) -> None:
    for cls, stats in by_class.items():
        lo, hi = bootstrap_recall_ci(stats["n_proteins"], stats["recall"])
        stats["recall_ci_lower"] = lo
        stats["recall_ci_upper"] = hi


def add_overall_ci(overall: dict) -> None:
    lo, hi = bootstrap_recall_ci(overall["n_proteins"], overall["recall"])
    overall["recall_ci_lower"] = lo
    overall["recall_ci_upper"] = hi


print("Processing binarized_fiveclass_eukaryote_eval_result.json ...")
with open(BINARIZED_JSON) as f:
    data = json.load(f)

add_class_cis(data["by_class"])

n_total   = data["n_algae_proteins"]
r_overall = data["overall_recall"]
lo, hi    = bootstrap_recall_ci(n_total, r_overall)
data["overall_recall_ci_lower"] = lo
data["overall_recall_ci_upper"] = hi

with open(BINARIZED_JSON, "w") as f:
    json.dump(data, f, indent=2)
print(f"  Saved to {BINARIZED_JSON}")

for cls, stats in data["by_class"].items():
    print(f"  {cls:<15} recall={stats['recall']:.4f}  "
          f"95% CI [{stats['recall_ci_lower']:.4f}, {stats['recall_ci_upper']:.4f}]")

print("\nDone.")
