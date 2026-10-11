# Homology-based inference vs. ProtT5 on the leakage-filtered algae: recall per class and across e-value thresholds.

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

BASE_DIR      = Path(__file__).resolve().parents[2]
HOMOLOGY_JSON = BASE_DIR / "results" / "homology_inference" / "homology_inference_result.json"
ML_JSON       = BASE_DIR / "results" / "training_results" / "binary_eukaryote_eval_clustered_result.json"
OUT_DIR       = BASE_DIR / "results" / "plots"
OUT_DIR.mkdir(parents=True, exist_ok=True)

CLASS_COLORS = {
    "thermophile":  "#E05C5C",
    "psychrophile": "#5C8AE0",
    "halophile":    "#E0A040",
    "acidophile":   "#6ABF69",
    "overall":      "#333333",
}

CLASS_ORDER = ["thermophile", "psychrophile", "halophile", "acidophile", "overall"]

plt.rcParams.update({
    "font.family":       "sans-serif",
    "font.size":         12,
    "axes.spines.top":   False,
    "axes.spines.right": False,
    "figure.dpi":        150,
})


def load_data():
    with open(HOMOLOGY_JSON) as f:
        homology = json.load(f)
    with open(ML_JSON) as f:
        ml = json.load(f)
    return homology, ml


def get_ml_recalls(ml):
    recalls = {}
    for cls in ["thermophile", "psychrophile", "halophile", "acidophile"]:
        if cls in ml["by_class"]:
            recalls[cls] = ml["by_class"][cls]["recall"]
    overall = ml.get("overall_recall") or ml.get("overall", {}).get("recall")
    if overall is not None:
        recalls["overall"] = overall
    return recalls


CLASS_LABELS = {
    "thermophile":  "Thermophile",
    "psychrophile": "Psychrophile",
    "halophile":    "Halophile",
    "acidophile":   "Acidophile",
    "overall":      "Overall",
}

THRESHOLD_LABELS = ["1e−1", "1e−3", "1e−5", "1e−10", "1e−20", "1e−50", "1e−100"]


def main():
    homology, ml = load_data()

    thresholds = [float(t) for t in homology["thresholds"]]
    results    = homology["results"]
    ml_recalls = get_ml_recalls(ml)

    first_result = list(results.values())[0]
    classes      = [c for c in CLASS_ORDER if c in first_result]

    fig, (ax_bar, ax_line) = plt.subplots(1, 2, figsize=(14, 6),
                                           gridspec_kw={"width_ratios": [1.1, 1]})

    # Left panel: homology at its best threshold vs. ProtT5
    best_thr    = thresholds[0]
    best_result = results[str(best_thr)]

    x      = np.arange(len(classes))
    width  = 0.35
    colors = [CLASS_COLORS.get(c, "#888") for c in classes]

    hom_vals = [best_result[c]["recall"] for c in classes]
    ml_vals  = [ml_recalls.get(c, 0.0)   for c in classes]

    bars_h = ax_bar.bar(x - width / 2, hom_vals, width,
                        color=colors, alpha=0.45, edgecolor="white",
                        label=f"Homology (e ≤ {best_thr})")
    bars_m = ax_bar.bar(x + width / 2, ml_vals, width,
                        color=colors, alpha=1.0, edgecolor="white",
                        label="ProtT5 FNN")

    for bar, val in zip(bars_h, hom_vals):
        ax_bar.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.015,
                    f"{val:.2f}", ha="center", va="bottom", fontsize=10, color="#555")
    for bar, val in zip(bars_m, ml_vals):
        ax_bar.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.015,
                    f"{val:.2f}", ha="center", va="bottom", fontsize=10, fontweight="bold")

    ax_bar.set_xticks(x)
    ax_bar.set_xticklabels([CLASS_LABELS[c] for c in classes], fontsize=11)
    ax_bar.set_ylim(0, 1.12)
    ax_bar.set_ylabel("Recall (fraction flagged as extremophile)", fontsize=11)
    ax_bar.set_title("Homology vs. ProtT5 — Best Threshold", fontsize=12)
    ax_bar.axhline(1.0, color="grey", linewidth=0.5, linestyle="--", alpha=0.4)
    ax_bar.legend(frameon=False, fontsize=10, loc="upper left")

    # Right panel: homology recall across e-value thresholds
    x_pos = np.arange(len(thresholds))

    for cls in classes:
        recalls = [results[str(t)][cls]["recall"] for t in thresholds]
        color   = CLASS_COLORS.get(cls, "#888888")
        lw      = 2.5 if cls == "overall" else 1.6
        ax_line.plot(x_pos, recalls, color=color, linewidth=lw,
                     marker="o", markersize=5, label=CLASS_LABELS[cls])

    ax_line.set_xticks(x_pos)
    ax_line.set_xticklabels(THRESHOLD_LABELS, fontsize=10, rotation=30, ha="right")
    ax_line.set_xlabel("E-value threshold  (stricter →)", fontsize=11)
    ax_line.set_ylabel("Recall", fontsize=11)
    ax_line.set_ylim(0, 0.30)
    ax_line.set_title("Homology Recall by Threshold", fontsize=12)
    ax_line.text(0.98, 0.97,
                 f"n = {homology['n_proteins']:,} proteins\n"
                 f"{homology['n_with_hits']:,} with ≥1 hit",
                 transform=ax_line.transAxes, ha="right", va="top",
                 fontsize=9, color="#888")

    fig.suptitle(
        "Homology-Based Inference vs. ProtT5 — Filtered Eukaryote Set",
        fontsize=13, y=1.01,
    )

    # Shared legend below both panels
    bar_handles = [
        mpatches.Patch(color="#888888", alpha=0.45, label=f"Homology (e ≤ 0.1)"),
        mpatches.Patch(color="#888888", alpha=1.00, label="ProtT5 FNN"),
    ]
    class_handles = [
        mpatches.Patch(color=CLASS_COLORS.get(c, "#888"), label=CLASS_LABELS[c])
        for c in classes
    ]
    fig.legend(handles=bar_handles + class_handles, frameon=False, fontsize=10,
               loc="lower center", ncol=len(bar_handles) + len(class_handles),
               bbox_to_anchor=(0.5, -0.06))

    plt.tight_layout()
    out = OUT_DIR / "homology_vs_ml_recall_clustered.png"
    plt.savefig(out, bbox_inches="tight")
    plt.close()
    print(f"Saved: {out}")


if __name__ == "__main__":
    main()
