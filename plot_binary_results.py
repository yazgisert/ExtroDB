"""
plot_binary_results.py
----------------------
Plots comparing binary ProtT5 classifier results against THERMO.

Run from thesis_data_collection/:
    python plot_binary_results.py
"""

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

OUT_DIR = Path("results/plots")
OUT_DIR.mkdir(parents=True, exist_ok=True)

plt.rcParams.update({
    "font.family":       "sans-serif",
    "font.size":         12,
    "axes.spines.top":   False,
    "axes.spines.right": False,
    "figure.dpi":        150,
})

# ── Load results ───────────────────────────────────────────────────────────────
SPLITS = ["split1_20pct", "split2_20pct", "split1_30pct", "split2_30pct"]
results = {}
for split in SPLITS:
    with open(f"results/training_results/{split}_binary_result.json") as f:
        results[split] = json.load(f)

# THERMO reference (from thermo_bacdive_result.json)
with open("results/training_results/thermo_bacdive_result.json") as f:
    thermo = json.load(f)

THERMO_METRICS = {
    "accuracy":           thermo["test_metrics"]["accuracy"],
    "macro-f1_score":     thermo["test_metrics"]["macro-f1_score"],
    "f1_thermophile":     thermo["test_metrics"]["f1_score"],      # thermophile F1
    "matthews-corr-coeff": thermo["test_metrics"]["matthews-corr-coeff"],
}

COLORS = {
    "split1_20pct": "#5C8AE0",
    "split2_20pct": "#A8C4F0",
    "split1_30pct": "#E0A85C",
    "split2_30pct": "#F0D4A8",
    "THERMO":       "#E05C5C",
}

LABELS = {
    "split1_20pct": "ProtT5 split1 20%",
    "split2_20pct": "ProtT5 split2 20%",
    "split1_30pct": "ProtT5 split1 30%",
    "split2_30pct": "ProtT5 split2 30%",
    "THERMO":       "THERMO (Pfam RF)",
}


# ── Plot 1: Side-by-side metric comparison ────────────────────────────────────
def plot_metric_comparison():
    metrics     = ["accuracy", "macro-f1_score", "f1_thermophile", "matthews-corr-coeff"]
    metric_labels = ["Accuracy", "Macro-F1", "F1 (Thermophile)", "MCC"]

    models = SPLITS + ["THERMO"]
    x      = np.arange(len(metrics))
    w      = 0.15
    offsets = np.linspace(-(len(models)-1)/2 * w, (len(models)-1)/2 * w, len(models))

    fig, ax = plt.subplots(figsize=(13, 5.5))

    for model, offset in zip(models, offsets):
        if model == "THERMO":
            vals = [THERMO_METRICS[m] for m in metrics]
        else:
            vals = [results[model]["test_metrics"][m] for m in metrics]
        ax.bar(x + offset, vals, w,
               label=LABELS[model], color=COLORS[model], edgecolor="white")

    ax.set_xticks(x)
    ax.set_xticklabels(metric_labels, fontsize=12)
    ax.set_ylim(0.5, 1.08)
    ax.set_ylabel("Score")
    ax.set_title("Binary Thermophile Classifier — ProtT5 vs. THERMO", fontsize=13)
    ax.legend(frameon=False, fontsize=9, loc="upper left",
              bbox_to_anchor=(0.0, 1.0), ncol=5)
    ax.axhline(1.0, color="grey", linewidth=0.5, linestyle="--", alpha=0.5)

    plt.tight_layout()
    out = OUT_DIR / "binary_metric_comparison.png"
    plt.savefig(out, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {out}")


# ── Plot 2: F1 thermophile vs non-thermophile per model ───────────────────────
def plot_f1_breakdown():
    models = SPLITS + ["THERMO"]
    f1_thermo    = []
    f1_nonthermo = []

    for model in models:
        if model == "THERMO":
            f1_thermo.append(THERMO_METRICS["f1_thermophile"])
            # THERMO non-thermo F1 derived from confusion matrix
            # TP=67 FN=4 FP=1 TN=166 → F1_nonthermo = 2*166/(2*166+4+1)
            f1_nonthermo.append(round(2*166/(2*166+4+1), 4))
        else:
            m = results[model]["test_metrics"]
            f1_thermo.append(m["f1_thermophile"])
            f1_nonthermo.append(m["f1_non_thermophile"])

    x  = np.arange(len(models))
    w  = 0.32
    fig, ax = plt.subplots(figsize=(11, 5.5))

    b1 = ax.bar(x - w/2, f1_thermo,    w, label="F1 — Thermophile",     color="#E05C5C", edgecolor="white")
    b2 = ax.bar(x + w/2, f1_nonthermo, w, label="F1 — Non-Thermophile", color="#5C8AE0", edgecolor="white")

    for bar in list(b1) + list(b2):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.008,
                f"{bar.get_height():.3f}", ha="center", va="bottom", fontsize=9)

    ax.set_xticks(x)
    ax.set_xticklabels([LABELS[m] for m in models], fontsize=10, rotation=12, ha="right")
    ax.set_ylim(0.5, 1.12)
    ax.set_ylabel("F1 Score")
    ax.set_title("Binary Classifier — F1 per Class", fontsize=13)
    ax.legend(frameon=False, loc="upper left")

    plt.tight_layout()
    out = OUT_DIR / "binary_f1_breakdown.png"
    plt.savefig(out, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {out}")


# ── Run all ────────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    print("Generating binary classifier plots ...")
    plot_metric_comparison()
    plot_f1_breakdown()
    print(f"\nDone. All plots saved to {OUT_DIR}")
