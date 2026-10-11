# Archaea bias plot from binary_domain_bias_result.json: MCC on all test proteins vs. bacteria only per embedder,
# and ProtT5 extremophile recall per class for archaeal vs. bacterial proteins.

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

RESULT = Path("results/training_results/binary_domain_bias_result.json")
OUT    = Path("results/plots/binary_domain_bias_clustered.png")

ALL_COLOR      = "#1B4F72"
BACTERIA_COLOR = "#7FB3D5"
ARCHAEA_COLOR  = "#C0392B"
EMBEDDERS = {"prott5": "ProtT5", "esm2_650m": "ESM2 650M", "esm2_8m": "ESM2 8M", "onehot": "One-hot"}
CLASSES   = ["acidophile", "alkaliphile", "halophile", "psychrophile", "thermophile"]

plt.rcParams.update({
    "font.family":       "sans-serif",
    "font.size":         12,
    "axes.spines.top":   False,
    "axes.spines.right": False,
    "figure.dpi":        150,
})


def bars_with_ci(ax, x, entries, color, label):
    values = np.array([e["value"] for e in entries])
    errors = np.array([values - [e["ci_lower"] for e in entries], [e["ci_upper"] for e in entries] - values])
    ax.bar(x, values, 0.38, color=color, edgecolor="white", label=label, yerr=errors, capsize=3,
           error_kw={"elinewidth": 1.1, "ecolor": "#333"})
    for xi, v, top in zip(x, values, values + errors[1]):
        ax.text(xi, top + 0.012, f"{v:.2f}", ha="center", va="bottom", fontsize=9)


if __name__ == "__main__":
    r = json.load(open(RESULT))
    fig, (left, right) = plt.subplots(1, 2, figsize=(13, 5), gridspec_kw={"width_ratios": [1.2, 1.4]})

    x = np.arange(len(EMBEDDERS))
    bars_with_ci(left, x - 0.2, [r["models"][e]["all"]["mcc"] for e in EMBEDDERS], ALL_COLOR, "All test proteins")
    bars_with_ci(left, x + 0.2, [r["models"][e]["bacteria_only"]["mcc"] for e in EMBEDDERS], BACTERIA_COLOR,
                 f"Bacteria only ({r['n_test_archaeal_extremophile']} archaeal extremophiles removed)")
    left.set_xticks(x)
    left.set_xticklabels(EMBEDDERS.values(), fontsize=10)
    left.set_ylim(0, 0.75)
    left.set_ylabel("MCC")
    left.set_title("Binary test set (clustered split)", fontsize=11)
    left.legend(fontsize=9, frameon=False, loc="upper right")

    per_class = r["prott5_recall_per_class"]
    x = np.arange(len(CLASSES))
    for offset, domain, color, label in [(-0.2, "archaea", ARCHAEA_COLOR, "Archaeal proteins"),
                                         (0.2, "bacteria", BACTERIA_COLOR, "Bacterial proteins")]:
        for xi, cls in zip(x, CLASSES):
            entry = per_class[cls].get(domain)
            if entry is None:
                right.text(xi + offset, 0.02, "none", ha="center", va="bottom", fontsize=8, color="#888", rotation=90)
                continue
            right.bar(xi + offset, entry["recall"], 0.38, color=color, edgecolor="white",
                      label=label if cls == "halophile" else None)
            right.text(xi + offset, entry["recall"] + 0.012, f"{entry['recall']:.2f}\nn={entry['n']}", ha="center",
                       va="bottom", fontsize=7.5)
    right.set_xticks(x)
    right.set_xticklabels([c.capitalize() for c in CLASSES], fontsize=10)
    right.set_ylim(0, 1.18)
    right.set_ylabel("Extremophile recall")
    right.set_title("ProtT5: recall per class and domain", fontsize=11)
    right.legend(fontsize=9, frameon=False, loc="upper left", ncol=2)

    fig.tight_layout()
    fig.savefig(OUT, bbox_inches="tight")
    print(f"Saved to {OUT}")
