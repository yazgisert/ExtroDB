# Plot of the transmembrane protein fraction per class, from tm_summary.json.

import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "plotting"))
from style import CLASS_COLORS as CLASS_COLORS_STYLE, NON_EXTREMOPHILE

SUMMARY_JSON = Path("results/tm_analysis/tm_summary.json")
OUT_DIR      = Path("results/plots")
OUT_DIR.mkdir(parents=True, exist_ok=True)

CLASS_COLORS = {
    "non_extremophile": NON_EXTREMOPHILE,
    "acidophile":        CLASS_COLORS_STYLE["acidophile"],
    "alkaliphile":       CLASS_COLORS_STYLE["alkaliphile"],
    "halophile":         CLASS_COLORS_STYLE["halophile"],
    "psychrophile":      CLASS_COLORS_STYLE["psychrophile"],
    "thermophile":       CLASS_COLORS_STYLE["thermophile"],
}

CLASS_LABELS = {
    "non_extremophile": "Non-\nextremophile",
    "acidophile":        "Acidophile",
    "alkaliphile":       "Alkaliphile",
    "halophile":         "Halophile",
    "psychrophile":      "Psychrophile",
    "thermophile":       "Thermophile",
}

plt.rcParams.update({
    "font.family":       "sans-serif",
    "font.size":         12,
    "axes.spines.top":   False,
    "axes.spines.right": False,
    "figure.dpi":        150,
})


def main():
    with open(SUMMARY_JSON) as f:
        data = json.load(f)

    # Build ordered list: non-extremophile first, then classes sorted by TM fraction
    ne = data["non_extremophile"]
    classes_sorted = sorted(
        data["extremophile_classes"].items(),
        key=lambda kv: kv[1]["tm_fraction_of_mapped"]
    )

    groups = [("non_extremophile", ne["tm_fraction"], ne["n_queried"], ne["n_with_tm"])]
    for cls, stats in classes_sorted:
        groups.append((
            cls,
            stats["tm_fraction_of_mapped"],
            stats["n_mapped_to_uniprot"],
            stats["n_with_tm"],
        ))

    labels  = [CLASS_LABELS[g[0]] for g in groups]
    fracs   = [g[1] * 100 for g in groups]          # convert to percent
    ns      = [g[2] for g in groups]
    n_tms   = [g[3] for g in groups]
    colors  = [CLASS_COLORS[g[0]] for g in groups]

    x = np.arange(len(groups))

    fig, ax = plt.subplots(figsize=(11, 6.0))
    fig.subplots_adjust(bottom=0.22, top=0.88)

    bars = ax.bar(x, fracs, color=colors, edgecolor="white", width=0.55)

    baseline = fracs[0]

    for bar, frac in zip(bars, fracs):
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height() + 0.6,
            f"{frac:.1f}%",
            ha="center", va="bottom", fontsize=10, fontweight="bold",
        )

    ax.axhline(baseline, color=NON_EXTREMOPHILE, linewidth=1.2, linestyle="--", alpha=0.6)
    ax.text(len(groups) - 0.52, baseline + 0.7,
            "non-extrem. baseline", fontsize=8, color=NON_EXTREMOPHILE,
            va="bottom", ha="right")

    combined_labels = [f"{lbl}\nn={n:,}" for lbl, n in zip(labels, ns)]
    ax.set_xticks(x)
    ax.set_xticklabels(combined_labels, fontsize=10, linespacing=1.6)

    ax.set_ylim(0, max(fracs) * 1.28)
    ax.set_ylabel("Transmembrane protein fraction (%)", fontsize=12)
    ax.set_title(
        "Transmembrane Protein Fraction by Extremophile Class",
        fontsize=14, pad=12,
    )

    # Methodology note as figure-level text, well below the axes
    fig.text(
        0.5, 0.01,
        "UniProt TM features  ·  extremophiles via WP_→UniProt mapping  ·  ~40% of WP_ IDs mapped",
        ha="center", va="bottom", fontsize=8, color="#aaa",
    )

    ax.spines["bottom"].set_visible(False)
    ax.tick_params(axis="x", length=0)

    plt.tight_layout()
    out = OUT_DIR / "tm_fraction_comparison.png"
    plt.savefig(out, bbox_inches="tight")
    plt.close()
    print(f"Saved: {out}")


if __name__ == "__main__":
    main()
