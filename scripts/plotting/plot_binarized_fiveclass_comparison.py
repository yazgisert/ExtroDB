# Plot of algae recall per class: binarized 5-class vs. dedicated binary model.

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

OUT_DIR       = Path("results/plots")
BINARY_JSON   = Path("results/training_results/binary_eukaryote_eval_clustered_result.json")
FIVECLASS_JSON = Path("results/training_results/binarized_fiveclass_eukaryote_eval_result.json")

OUT_DIR.mkdir(parents=True, exist_ok=True)

CLASS_COLORS = {
    "thermophile":  "#E05C5C",
    "psychrophile": "#5C8AE0",
    "halophile":    "#E0A040",
    "acidophile":   "#6ABF69",
}

CLASS_ORDER = ["halophile", "psychrophile", "acidophile", "thermophile"]

plt.rcParams.update({
    "font.family":       "sans-serif",
    "font.size":         12,
    "axes.spines.top":   False,
    "axes.spines.right": False,
    "figure.dpi":        150,
})


def main():
    with open(BINARY_JSON) as f:
        binary = json.load(f)
    with open(FIVECLASS_JSON) as f:
        fiveclass = json.load(f)

    binary_by_class    = binary["by_class"]
    fiveclass_by_class = fiveclass["by_class"]

    binary_recalls    = [binary_by_class[c]["recall"]    for c in CLASS_ORDER]
    fiveclass_recalls = [fiveclass_by_class[c]["recall"] for c in CLASS_ORDER]

    def get_yerr(by_class, recalls):
        lo_errs = [r - by_class[c].get("recall_ci_lower", r) for c, r in zip(CLASS_ORDER, recalls)]
        hi_errs = [by_class[c].get("recall_ci_upper", r) - r for c, r in zip(CLASS_ORDER, recalls)]
        return np.array([lo_errs, hi_errs])

    binary_yerr    = get_yerr(binary_by_class, binary_recalls)
    fiveclass_yerr = get_yerr(fiveclass_by_class, fiveclass_recalls)

    x     = np.arange(len(CLASS_ORDER))
    width = 0.35

    fig, ax = plt.subplots(figsize=(9, 5))

    bars_binary = ax.bar(x - width / 2, binary_recalls, width,
                         color=[CLASS_COLORS[c] for c in CLASS_ORDER],
                         alpha=1.0, edgecolor="white",
                         yerr=binary_yerr, capsize=4,
                         error_kw={"elinewidth": 1.2, "ecolor": "#222", "capthick": 1.2})

    bars_five = ax.bar(x + width / 2, fiveclass_recalls, width,
                       color=[CLASS_COLORS[c] for c in CLASS_ORDER],
                       alpha=0.35, edgecolor="white",
                       yerr=fiveclass_yerr, capsize=4,
                       error_kw={"elinewidth": 1.2, "ecolor": "#555", "capthick": 1.2})

    for bar, val in zip(bars_binary, binary_recalls):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.025,
                f"{val:.3f}", ha="center", va="bottom", fontsize=10, fontweight="bold")

    for bar, val in zip(bars_five, fiveclass_recalls):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.025,
                f"{val:.3f}", ha="center", va="bottom", fontsize=10, color="#555")

    ax.set_xticks(x)
    ax.set_xticklabels([c.capitalize() for c in CLASS_ORDER], fontsize=12)
    ax.set_ylim(0, 1.15)
    ax.set_ylabel("Recall (fraction flagged as extremophile)")
    ax.set_title(
        "Eukaryote Recall — Dedicated Binary Classifier vs. Binarized 5-Class Model",
        fontsize=12
    )
    ax.axhline(1.0, color="grey", linewidth=0.5, linestyle="--", alpha=0.4)

    legend_handles = [
        mpatches.Patch(color="#888888", alpha=1.0,  label="Dedicated binary classifier"),
        mpatches.Patch(color="#888888", alpha=0.35, label="Binarized 5-class (max prob ≥ 0.5)"),
    ]
    ax.legend(handles=legend_handles, frameon=False, fontsize=9.5, loc="lower right")

    plt.tight_layout()
    out = OUT_DIR / "binarized_fiveclass_vs_binary_clustered.png"
    plt.savefig(out, bbox_inches="tight")
    plt.close()
    print(f"Saved: {out}")


if __name__ == "__main__":
    main()
