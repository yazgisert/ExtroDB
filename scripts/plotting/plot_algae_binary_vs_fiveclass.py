# Plot of the binary model's algae recall vs. the 5-class model's algae results per class.

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

BINARY_JSON    = Path("results/training_results/binary_eukaryote_eval_clustered_result.json")
FIVECLASS_JSON = Path("results/training_results/eukaryote_fiveclass_filtered_result.json")
OUT            = Path("results/plots/algae_binary_vs_fiveclass_clustered.png")

CLASS_COLORS = {
    "thermophile":  "#E05C5C",
    "psychrophile": "#5C8AE0",
    "halophile":    "#E0A040",
    "acidophile":   "#6ABF69",
}

plt.rcParams.update({
    "font.family":       "sans-serif",
    "font.size":         12,
    "axes.spines.top":   False,
    "axes.spines.right": False,
    "figure.dpi":        150,
})


if __name__ == "__main__":
    binary = json.load(open(BINARY_JSON))["by_class"]
    fiveclass = json.load(open(FIVECLASS_JSON))["eukaryote_eval"]
    classes = list(CLASS_COLORS)

    recall = np.array([binary[c]["recall"] for c in classes])
    recall_err = np.array([[binary[c]["recall"] - binary[c]["recall_ci_lower"] for c in classes],
                           [binary[c]["recall_ci_upper"] - binary[c]["recall"] for c in classes]])
    f1 = np.array([fiveclass["per_class_f1"][c]["value"] for c in classes])
    f1_err = np.array([[fiveclass["per_class_f1"][c]["value"] - fiveclass["per_class_f1"][c]["ci_lower"] for c in classes],
                       [fiveclass["per_class_f1"][c]["ci_upper"] - fiveclass["per_class_f1"][c]["value"] for c in classes]])

    x, width = np.arange(len(classes)), 0.35
    colors = [CLASS_COLORS[c] for c in classes]
    fig, ax = plt.subplots(figsize=(9, 5))
    ax.bar(x - width / 2, recall, width, color=colors, edgecolor="white", yerr=recall_err, capsize=4,
           error_kw={"elinewidth": 1.2, "ecolor": "#222"})
    ax.bar(x + width / 2, f1, width, color=colors, alpha=0.45, edgecolor="white", yerr=f1_err, capsize=4,
           error_kw={"elinewidth": 1.2, "ecolor": "#555"})
    for xi, r, f in zip(x, recall, f1):
        ax.text(xi - width / 2, r + 0.025, f"{r:.2f}", ha="center", fontweight="bold")
        ax.text(xi + width / 2, f + 0.025, f"{f:.2f}", ha="center", fontsize=11)

    ax.set_xticks(x)
    ax.set_xticklabels([c.capitalize() for c in classes])
    ax.set_ylim(0, 1.15)
    ax.set_ylabel("Score")
    ax.axhline(1.0, color="grey", linewidth=0.5, linestyle="--", alpha=0.4)
    ax.bar(0, 0, color="#555", label="Binary recall (extremophile vs. non-extremophile)")
    ax.bar(0, 0, color="#555", alpha=0.45, label="5-class F1 (which type of extremophile)")
    ax.legend(frameon=False, fontsize=10, loc="upper right")
    ax.set_title(f"Leakage-Filtered Algae (n = {fiveclass['n']:,}) — Binary Recall vs. 5-Class F1 (ProtT5)", fontsize=12)

    plt.tight_layout()
    plt.savefig(OUT, bbox_inches="tight")
    plt.close()
    print(f"Saved: {OUT}")
