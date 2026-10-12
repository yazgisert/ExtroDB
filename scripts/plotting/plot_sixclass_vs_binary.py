# Plot of the binarized 6-class model vs. the dedicated binary model on the clustered test set.

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from style import MAIN, SECONDARY

OUT_DIR       = Path("results/plots")
SIXCLASS_JSON = Path("results/training_results/sixclass_clustered_result.json")

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

BINARY_COLOR   = MAIN
SIXCLASS_COLOR = SECONDARY

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
        result = json.load(f)

    binary_vals, binary_err, binary_top = values_and_errors(result["dedicated_binary"])
    six_vals, six_err, six_top          = values_and_errors(result["binarized_sixclass"])

    x     = np.arange(len(METRICS))
    width = 0.38

    fig, ax = plt.subplots(figsize=(11, 5))

    ax.bar(x - width / 2, binary_vals, width, color=BINARY_COLOR, edgecolor="white",
           yerr=binary_err, capsize=4, error_kw={"elinewidth": 1.2, "ecolor": "#222"},
           label="Dedicated binary classifier")
    ax.bar(x + width / 2, six_vals, width, color=SIXCLASS_COLOR, edgecolor="white",
           yerr=six_err, capsize=4, error_kw={"elinewidth": 1.2, "ecolor": "#555"},
           label="Binarized 6-class (P(non-extremophile) < 0.5)")

    for xi, val, top in zip(x - width / 2, binary_vals, binary_top):
        ax.text(xi, top + 0.01, f"{val:.3f}", ha="center", va="bottom", fontsize=9, fontweight="bold")
    for xi, val, top in zip(x + width / 2, six_vals, six_top):
        ax.text(xi, top + 0.01, f"{val:.3f}", ha="center", va="bottom", fontsize=9, color="#444")

    ax.set_xticks(x)
    ax.set_xticklabels(METRICS.values(), fontsize=11)
    ax.set_ylim(0, 1.1)
    ax.set_ylabel("Score on prokaryote test set")
    ax.set_title(f"Extremophile Detection — Dedicated Binary vs. Binarized 6-Class Model "
                 f"(n = {result['n_test']:,})", fontsize=12)
    ax.legend(frameon=False, fontsize=10, loc="upper center",
              bbox_to_anchor=(0.5, -0.15), ncol=2)

    plt.tight_layout()
    out = OUT_DIR / "sixclass_vs_binary_clustered.png"
    plt.savefig(out, bbox_inches="tight")
    plt.close()
    print(f"Saved: {out}")


if __name__ == "__main__":
    main()
