# Plot of binary algae recall before vs. after the leakage filter.

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

OUT_DIR      = Path("results/plots")
ORIGINAL_JSON = Path("results/training_results/binary_eukaryote_eval_unfiltered_clustered_result.json")
FILTERED_JSON = Path("results/training_results/binary_eukaryote_eval_clustered_result.json")

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
    with open(ORIGINAL_JSON) as f:
        original = json.load(f)
    with open(FILTERED_JSON) as f:
        filtered = json.load(f)

    orig_by_class = original["by_class"]
    filt_by_class = filtered["by_class"]

    classes      = CLASS_ORDER
    orig_recalls = [orig_by_class[c]["recall"] for c in classes]
    filt_recalls = [filt_by_class[c]["recall"] for c in classes]
    filt_ns      = [filt_by_class[c]["n_proteins"] for c in classes]
    orig_ns      = [orig_by_class[c]["n_proteins"] for c in classes]

    def get_yerr(by_class, recalls):
        lo_errs = [r - by_class[c].get("recall_ci_lower", r) for c, r in zip(classes, recalls)]
        hi_errs = [by_class[c].get("recall_ci_upper", r) - r for c, r in zip(classes, recalls)]
        return np.array([lo_errs, hi_errs])

    orig_yerr = get_yerr(orig_by_class, orig_recalls)
    filt_yerr = get_yerr(filt_by_class, filt_recalls)

    x     = np.arange(len(classes))
    width = 0.35

    fig, ax = plt.subplots(figsize=(9, 5))

    bars_orig = ax.bar(x - width / 2, orig_recalls, width,
                       color=[CLASS_COLORS[c] for c in classes],
                       alpha=0.40, edgecolor="white",
                       yerr=orig_yerr, capsize=4,
                       error_kw={"elinewidth": 1.2, "ecolor": "#555", "capthick": 1.2})

    bars_filt = ax.bar(x + width / 2, filt_recalls, width,
                       color=[CLASS_COLORS[c] for c in classes],
                       alpha=1.00, edgecolor="white",
                       yerr=filt_yerr, capsize=4,
                       error_kw={"elinewidth": 1.2, "ecolor": "#222", "capthick": 1.2})

    for bar, val, n_orig in zip(bars_orig, orig_recalls, orig_ns):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.025,
                f"{val:.3f}", ha="center", va="bottom", fontsize=10, color="#666")

    for bar, val, n_filt in zip(bars_filt, filt_recalls, filt_ns):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.025,
                f"{val:.3f}", ha="center", va="bottom", fontsize=10, fontweight="bold")

    ax.set_xticks(x)
    ax.set_xticklabels([c.capitalize() for c in classes], fontsize=12)
    ax.set_ylim(0, 1.15)
    ax.set_ylabel("Recall (fraction flagged as extremophile)")
    ax.set_title(
        "Binary Classifier on Eukaryote Algae — Before vs. After Leakage Filtering",
        fontsize=12
    )
    ax.axhline(1.0, color="grey", linewidth=0.5, linestyle="--", alpha=0.4)

    legend_handles = [
        mpatches.Patch(color="#888888", alpha=0.40,
                       label=f"Before filtering  (n = {original['n_algae_proteins']:,})"),
        mpatches.Patch(color="#888888", alpha=1.00,
                       label=f"After filtering   (n = {filtered['n_algae_proteins']:,}, removed proteins ≥ 30% identity to train+val)"),
    ]
    ax.legend(handles=legend_handles, frameon=False, fontsize=9.5, loc="lower right")

    plt.tight_layout()
    out = OUT_DIR / "binary_eukaryote_filtered_comparison_clustered.png"
    plt.savefig(out, bbox_inches="tight")
    plt.close()
    print(f"Saved: {out}")


if __name__ == "__main__":
    main()
