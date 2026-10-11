# Plot of homology-based inference vs. ProtT5 on the 5-class test set.

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

RESULT = Path("results/training_results/homology_fiveclass_result.json")
OUT    = Path("results/plots/homology_fiveclass_vs_prott5.png")

HOMOLOGY_COLOR = "#A8A8A8"
PROTT5_COLOR   = "#1B4F72"
CLASSES = ["acidophile", "alkaliphile", "halophile", "psychrophile", "thermophile"]

plt.rcParams.update({
    "font.family":       "sans-serif",
    "font.size":         12,
    "axes.spines.top":   False,
    "axes.spines.right": False,
    "figure.dpi":        150,
})


def bars_with_ci(ax, x, entries, color, label, bold=False):
    values = np.array([e["value"] for e in entries])
    errors = np.array([values - [e["ci_lower"] for e in entries], [e["ci_upper"] for e in entries] - values])
    ax.bar(x, values, 0.38, color=color, edgecolor="white", label=label, yerr=errors, capsize=3,
           error_kw={"elinewidth": 1.1, "ecolor": "#333"})
    for xi, v, top in zip(x, values, values + errors[1]):
        ax.text(xi, top + 0.015, f"{v:.2f}", ha="center", va="bottom", fontsize=9,
                fontweight="bold" if bold else "normal")


if __name__ == "__main__":
    r = json.load(open(RESULT))
    fig, (left, right) = plt.subplots(1, 2, figsize=(13, 5), gridspec_kw={"width_ratios": [3, 1.3]})

    groups = ["Accuracy"] + [c.capitalize() for c in CLASSES]
    x = np.arange(len(groups))
    for offset, key, color, label, bold in [(-0.2, "homology", HOMOLOGY_COLOR, "Homology-based inference (best MMseqs2 hit)", False),
                                            (0.2, "prott5", PROTT5_COLOR, "ProtT5 + FNN", True)]:
        entries = [r[key]["accuracy"]] + [r[key]["per_class_recall"][c] for c in CLASSES]
        bars_with_ci(left, x + offset, entries, color, label, bold)
    left.axvline(0.5, color="#ccc", linewidth=1)
    left.set_xticks(x)
    left.set_xticklabels(groups, fontsize=10)
    left.set_ylim(0, 1.05)
    left.set_ylabel("Accuracy / per-class recall")
    left.set_title("All 5-class test proteins", fontsize=11)

    s = r["subsets"]
    labels = [f"With a hit\n(n = {s['with_hit']['n']:,})", f"Without a hit\n(n = {s['without_hit']['n']:,})"]
    right.bar([-0.2], [s["with_hit"]["homology_accuracy"]], 0.38, color=HOMOLOGY_COLOR)
    right.bar([0.2, 1.0], [s["with_hit"]["prott5_accuracy"], s["without_hit"]["prott5_accuracy"]], 0.38, color=PROTT5_COLOR)
    for xi, v, bold in [(-0.2, s["with_hit"]["homology_accuracy"], False), (0.2, s["with_hit"]["prott5_accuracy"], True),
                        (1.0, s["without_hit"]["prott5_accuracy"], True)]:
        right.text(xi, v + 0.015, f"{v:.2f}", ha="center", va="bottom", fontsize=9, fontweight="bold" if bold else "normal")
    right.text(0.8, 0.03, "no\nprediction", ha="center", va="bottom", fontsize=8, color="#888")
    right.set_xticks([0, 0.9])
    right.set_xticklabels(labels, fontsize=10)
    right.set_ylim(0, 1.05)
    right.set_ylabel("Accuracy")
    right.set_title("Split by homology to training", fontsize=11)

    handles, names = left.get_legend_handles_labels()
    fig.legend(handles, names, loc="lower center", ncol=2, frameon=False, fontsize=10, bbox_to_anchor=(0.5, -0.06))
    fig.suptitle(f"Homology-Based Inference vs. ProtT5 — Prokaryote 5-Class Test Set (n = {r['n_test']:,}, split1 30%)",
                 fontsize=13)
    plt.tight_layout()
    plt.savefig(OUT, bbox_inches="tight")
    plt.close()
    print(f"Saved: {OUT}")
