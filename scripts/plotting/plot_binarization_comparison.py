# Plot of the binarized 5-class model, the binarized 6-class model and the dedicated binary model on the clustered
# test set.

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from style import MAIN, RAMP, SECONDARY

OUT_DIR        = Path("results/plots")
SIXCLASS_JSON  = Path("results/training_results/sixclass_clustered_result.json")
FIVECLASS_JSON = Path("results/training_results/fiveclass_binarized_clustered_result.json")

OUT_DIR.mkdir(parents=True, exist_ok=True)

METRICS = {
    "accuracy":          "Accuracy",
    "balanced_accuracy": "Balanced\naccuracy",
    "precision":         "Precision",
    "recall":            "Recall",
    "specificity":       "Specificity",
    "f1":                "F1",
    "mcc":               "MCC",
}

plt.rcParams.update({
    "font.family":       "sans-serif",
    "font.size":         12,
    "axes.spines.top":   False,
    "axes.spines.right": False,
    "figure.dpi":        150,
})


def values_and_errors(results):
    values = np.array([results[m]["value"] for m in METRICS])
    lower  = np.array([results[m]["ci_lower"] for m in METRICS])
    upper  = np.array([results[m]["ci_upper"] for m in METRICS])
    return values, np.array([values - lower, upper - values]), upper


def main():
    with open(SIXCLASS_JSON) as f:
        sixclass = json.load(f)
    with open(FIVECLASS_JSON) as f:
        fiveclass = json.load(f)

    models = [
        ("Binarized 5-class (max prob ≥ 0.5)",          fiveclass["binarized_fiveclass"], RAMP[0], "#555"),
        ("Binarized 6-class (P(non-extremophile) < 0.5)", sixclass["binarized_sixclass"],  SECONDARY, "#555"),
        ("Dedicated binary classifier",                 sixclass["dedicated_binary"],     MAIN, "#222"),
    ]

    x     = np.arange(len(METRICS))
    width = 0.27

    fig, ax = plt.subplots(figsize=(12, 5))

    for i, (label, results, color, err_color) in enumerate(models):
        values, errors, top = values_and_errors(results)
        offset = (i - 1) * width
        ax.bar(x + offset, values, width, color=color, edgecolor="white",
               yerr=errors, capsize=3, error_kw={"elinewidth": 1.1, "ecolor": err_color},
               label=label)
        for xi, val, t in zip(x + offset, values, top):
            ax.text(xi, t + 0.01, f"{val:.2f}", ha="center", va="bottom", fontsize=8.5,
                    fontweight="bold" if i == 2 else "normal", color="#222" if i == 2 else "#444")

    ax.set_xticks(x)
    ax.set_xticklabels(METRICS.values(), fontsize=11)
    ax.set_ylim(0, 1.1)
    ax.set_ylabel("Score on prokaryote test set")
    ax.set_title(f"Extremophile Detection — 5-Class vs. 6-Class vs. Binary Model "
                 f"(same test set, n = {sixclass['n_test']:,})", fontsize=12)
    ax.legend(frameon=False, fontsize=10, loc="upper center",
              bbox_to_anchor=(0.5, -0.15), ncol=3)

    plt.tight_layout()
    out = OUT_DIR / "binarization_comparison_clustered.png"
    plt.savefig(out, bbox_inches="tight")
    plt.close()
    print(f"Saved: {out}")


if __name__ == "__main__":
    main()
