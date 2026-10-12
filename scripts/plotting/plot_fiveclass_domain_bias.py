# Archaea check plot for the 5-class ProtT5 model, from fiveclass_domain_bias_result.json: recall per class for
# archaeal vs. bacterial test proteins, and halophile / thermophile recall per subtype and domain.

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

RESULT = Path("results/training_results/fiveclass_domain_bias_result.json")
OUT    = Path("results/plots/fiveclass_domain_bias.png")

ARCHAEA_COLOR  = "#C0392B"
BACTERIA_COLOR = "#7FB3D5"
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
            if entry["recall"] is None:
                left.text(xi + offset, 0.02, "none", ha="center", va="bottom", fontsize=8, color="#888", rotation=90)
                continue
            left.bar(xi + offset, entry["recall"], 0.38, color=color, edgecolor="white",
                     label=label if cls == "halophile" else None)
            left.text(xi + offset, entry["recall"] + 0.012, f"{entry['recall']:.2f}\nn={entry['n']:,}", ha="center",
                      va="bottom", fontsize=7.5)
    left.set_xticks(x)
    left.set_xticklabels([c.capitalize() for c in CLASSES], fontsize=10)
    left.set_ylim(0, 1.18)
    left.set_ylabel("Recall")
    left.set_title(f"Recall per class (macro-F1 {r['all']['macro_f1']['value']:.3f} all, "
                   f"{r['bacteria_only']['macro_f1']['value']:.3f} bacteria only)", fontsize=11)
    left.legend(fontsize=9, frameon=False, loc="upper left", ncol=2)

    entries = [r["recall_by_subtype_and_domain"][key] for key, _ in SUBTYPES]
    colors = [ARCHAEA_COLOR if key.endswith("archaea") else BACTERIA_COLOR for key, _ in SUBTYPES]
    x = np.arange(len(SUBTYPES))
    right.bar(x, [e["recall"] for e in entries], 0.6, color=colors, edgecolor="white")
    for xi, e in zip(x, entries):
        right.text(xi, e["recall"] + 0.012, f"{e['recall']:.2f}\nn={e['n']:,}", ha="center", va="bottom", fontsize=7.5)
    right.axvline(3.5, color="#999", linewidth=0.8, linestyle="--")
    right.set_xticks(x)
    right.set_xticklabels([label for _, label in SUBTYPES], fontsize=8.5)
    right.set_ylim(0, 1.18)
    right.set_ylabel("Recall of the correct class")
    right.set_title("Halophiles and thermophiles by subtype and domain", fontsize=11)

    fig.tight_layout()
    fig.savefig(OUT, bbox_inches="tight")
    print(f"Saved to {OUT}")
