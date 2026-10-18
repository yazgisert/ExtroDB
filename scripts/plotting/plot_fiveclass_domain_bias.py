# Archaea check plots for the 5-class model: accuracy per class for archaeal vs. bacterial test proteins and per
# halophile / thermophile subtype (ProtT5), and the subtype accuracy of ProtT5 next to amino acid composition.

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from style import EMBEDDER_COLORS, MAIN, SECONDARY

RESULT        = Path("results/training_results/fiveclass_domain_bias_result.json")
RESULT_ONEHOT = Path("results/training_results/fiveclass_domain_bias_onehot_result.json")
OUT           = Path("results/plots/fiveclass_domain_bias.png")
OUT_EMBEDDERS = Path("results/plots/fiveclass_domain_bias_embedders.png")

ARCHAEA_COLOR  = MAIN
BACTERIA_COLOR = SECONDARY
CLASSES = ["acidophile", "alkaliphile", "halophile", "psychrophile", "thermophile"]
SUBTYPES = [
    ("extreme_halophile|archaea",  "Extreme\nhalophile\narchaea"),
    ("moderate_halophile|archaea", "Moderate\nhalophile\narchaea"),
    ("moderate_halophile|bacteria", "Moderate\nhalophile\nbacteria"),
    ("slight_halophile|bacteria",   "Slight\nhalophile\nbacteria"),
    ("hyperthermophile|archaea",    "Hyper-\nthermophile\narchaea"),
    ("thermophile|archaea",         "Thermophile\n\narchaea"),
    ("thermophile|bacteria",        "Thermophile\n\nbacteria"),
]

plt.rcParams.update({
    "font.family":       "sans-serif",
    "font.size":         12,
    "axes.spines.top":   False,
    "axes.spines.right": False,
    "figure.dpi":        150,
})


if __name__ == "__main__":
    r = json.load(open(RESULT))
    fig, (left, right) = plt.subplots(1, 2, figsize=(14, 5), gridspec_kw={"width_ratios": [1.1, 1.5]})

    x = np.arange(len(CLASSES))
    for offset, domain, color, label in [(-0.2, "archaea", ARCHAEA_COLOR, "Archaeal proteins"),
                                         (0.2, "bacteria", BACTERIA_COLOR, "Bacterial proteins")]:
        for xi, cls in zip(x, CLASSES):
            entry = r["per_class"][cls][domain]
            if entry["accuracy"] is None:
                left.text(xi + offset, 0.02, "none", ha="center", va="bottom", fontsize=8, color="#888", rotation=90)
                continue
            left.bar(xi + offset, entry["accuracy"], 0.38, color=color, edgecolor="white",
                     label=label if cls == "halophile" else None)
            left.text(xi + offset, entry["accuracy"] + 0.012, f"{entry['accuracy']:.2f}\nn={entry['n']:,}", ha="center",
                      va="bottom", fontsize=7.5)
    left.set_xticks(x)
    left.set_xticklabels([c.capitalize() for c in CLASSES], fontsize=10)
    left.set_ylim(0, 1.18)
    left.set_ylabel("Accuracy")
    left.set_title(f"Accuracy per class (overall: archaea {r['archaea_only']['accuracy']['value']:.2f}, "
                   f"bacteria {r['bacteria_only']['accuracy']['value']:.2f})", fontsize=11)
    left.legend(fontsize=9, frameon=False, loc="upper left", ncol=2)

    entries = [r["accuracy_by_subtype_and_domain"][key] for key, _ in SUBTYPES]
    colors = [ARCHAEA_COLOR if key.endswith("archaea") else BACTERIA_COLOR for key, _ in SUBTYPES]
    x = np.arange(len(SUBTYPES))
    right.bar(x, [e["accuracy"] for e in entries], 0.6, color=colors, edgecolor="white")
    for xi, e in zip(x, entries):
        right.text(xi, e["accuracy"] + 0.012, f"{e['accuracy']:.2f}\nn={e['n']:,}", ha="center", va="bottom", fontsize=7.5)
    right.axvline(3.5, color="#999", linewidth=0.8, linestyle="--")
    right.set_xticks(x)
    right.set_xticklabels([label for _, label in SUBTYPES], fontsize=8.5)
    right.set_ylim(0, 1.18)
    right.set_ylabel("Accuracy")
    right.set_title("Halophiles and thermophiles by subtype and domain", fontsize=11)

    fig.tight_layout()
    fig.savefig(OUT, bbox_inches="tight")
    print(f"Saved to {OUT}")

    # Same subtypes, ProtT5 next to amino acid composition
    onehot = json.load(open(RESULT_ONEHOT))
    fig, ax = plt.subplots(figsize=(11, 4.5))
    x = np.arange(len(SUBTYPES))
    for offset, res, color, label in [(-0.2, r, EMBEDDER_COLORS["prott5"], "ProtT5"),
                                      (0.2, onehot, EMBEDDER_COLORS["onehot"], "Composition")]:
        values = [res["accuracy_by_subtype_and_domain"][key]["accuracy"] for key, _ in SUBTYPES]
        ax.bar(x + offset, values, 0.38, color=color, edgecolor="white", label=label)
        for xi, v in zip(x, values):
            ax.text(xi + offset, v + 0.012, f"{v:.2f}", ha="center", va="bottom", fontsize=8)
    ax.axvline(3.5, color="#999", linewidth=0.8, linestyle="--")
    ax.set_xticks(x)
    ax.set_xticklabels([label for _, label in SUBTYPES], fontsize=8.5)
    ax.set_ylim(0, 1.1)
    ax.set_ylabel("Accuracy")
    ax.set_title("5-class model: halophiles and thermophiles by subtype and domain", fontsize=11)
    ax.legend(fontsize=9, frameon=False, loc="upper right", ncol=2)
    fig.tight_layout()
    fig.savefig(OUT_EMBEDDERS, bbox_inches="tight")
    print(f"Saved to {OUT_EMBEDDERS}")
