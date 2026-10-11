# THERMO on the 13 algae: predicted thermophile probability per species, coloured by true class.

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

OUT_DIR  = Path("results/plots")
RESULT   = Path("results/training_results/thermo_eukaryote_eval_result.json")

OUT_DIR.mkdir(parents=True, exist_ok=True)

# Palette consistent with the rest of the thesis plots
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


def shorten(name: str) -> str:
    return (name
            .replace("Nannochloropsis", "N.")
            .replace("Chlamydomonas sp.", "Chlamy.")
            .replace("Chlamydomonas", "Chlamy.")
            .replace("Fragilariopsis", "Frag.")
            .replace("Coccomyxa", "Coc.")
            .replace("Galdieria", "G.")
            .replace("Cyanidium", "Cya.")
            .replace("Dunaliella", "Dun.")
            .replace("Polarella", "Pol."))


def main():
    with open(RESULT) as f:
        data = json.load(f)

    records = data["per_proteome"]

    ordered = sorted(records, key=lambda r: r["proba_thermophile"], reverse=True)
    thermo  = [r for r in ordered if r["true_label"] == 1]
    other   = [r for r in ordered if r["true_label"] == 0]

    labels = [shorten(r["species"]) for r in ordered]
    probs  = [r["proba_thermophile"] for r in ordered]
    colors = [CLASS_COLORS[r["broad_class"]] for r in ordered]
    edge   = ["#333" if r["true_label"] == 1 else "white" for r in ordered]

    y = np.arange(len(ordered))

    fig, ax = plt.subplots(figsize=(9, 6))

    bars = ax.barh(y, probs, color=colors, edgecolor=edge,
                   linewidth=1.4, height=0.65)

    ax.axvline(0.5, color="#555", linewidth=1.2, linestyle="--",
               label="Decision threshold (0.5)")

    for i, (prob, r) in enumerate(zip(probs, ordered)):
        ax.text(prob + 0.006, i, f"{prob:.3f}",
                va="center", ha="left", fontsize=9,
                color="#444")

    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=10)
    ax.invert_yaxis()
    ax.set_xlim(0, 0.65)
    ax.set_xlabel("Predicted thermophile probability (THERMO model)")

    legend_patches = [mpatches.Patch(color=CLASS_COLORS[c], label=c.capitalize())
                      for c in ["thermophile", "psychrophile", "halophile", "acidophile"]]
    legend_patches.append(
        plt.Line2D([0], [0], color="#555", linewidth=1.2, linestyle="--",
                   label="Decision threshold (0.5)")
    )
    ax.legend(handles=legend_patches, frameon=False, fontsize=9,
              bbox_to_anchor=(1.01, 1), loc="upper left")

    fig.suptitle(
        "THERMO Model (Prokaryote-Trained) Applied to Eukaryote Algae\n"
        "Recall = 0/4 thermophile species  ·  Specificity = 9/9 non-thermophile species",
        fontsize=12, y=1.02
    )

    plt.tight_layout()
    out = OUT_DIR / "thermo_eukaryote_eval.png"
    plt.savefig(out, bbox_inches="tight")
    plt.close()
    print(f"Saved: {out}")


if __name__ == "__main__":
    main()
