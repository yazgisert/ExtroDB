# Binary thermophile classifier (four ProtT5 splits) vs. THERMO: metric comparison and F1 per class.

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from style import BASELINE, CLASS_COLORS, NON_EXTREMOPHILE, RAMP

OUT_DIR = Path("results/plots")
OUT_DIR.mkdir(parents=True, exist_ok=True)

plt.rcParams.update({
    "font.family":       "sans-serif",
    "font.size":         12,
    "axes.spines.top":   False,
    "axes.spines.right": False,
    "figure.dpi":        150,
})

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
    "split1_20pct": RAMP[1],
    "split2_20pct": RAMP[0],
    "split1_30pct": RAMP[3],
    "split2_30pct": RAMP[2],
    "THERMO":       BASELINE,
}

LABELS = {
    "split1_20pct": "ProtT5 split1 20%",
    "split2_20pct": "ProtT5 split2 20%",
    "split1_30pct": "ProtT5 split1 30%",
    "split2_30pct": "ProtT5 split2 30%",
    "THERMO":       "THERMO (Pfam RF)",
}


# Plot 1: metric comparison
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
            vals  = [THERMO_METRICS[m] for m in metrics]
            yerr  = None
        else:
            b     = results[model]["test_metrics_bootstrap"]
            vals  = [results[model]["test_metrics"][m] for m in metrics]
            lo    = [vals[i] - b[mk]["lower"] for i, mk in enumerate(metrics)]
            hi    = [b[mk]["upper"] - vals[i]  for i, mk in enumerate(metrics)]
            yerr  = [lo, hi]
        ax.bar(x + offset, vals, w,
               label=LABELS[model], color=COLORS[model], edgecolor="white",
               yerr=yerr, capsize=3, error_kw={"elinewidth": 1.0, "ecolor": "#555"})

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


# Plot 2: F1 per class and accuracy per model
def plot_f1_breakdown():
    models = SPLITS + ["THERMO"]
    f1_thermo    = []
    f1_nonthermo = []
    accuracies   = []

    for model in models:
        if model == "THERMO":
            f1_thermo.append(THERMO_METRICS["f1_thermophile"])
            # THERMO non-thermo F1 derived from confusion matrix
            # TP=67 FN=4 FP=1 TN=166 → F1_nonthermo = 2*166/(2*166+4+1)
            f1_nonthermo.append(round(2*166/(2*166+4+1), 4))
            accuracies.append(THERMO_METRICS["accuracy"])
        else:
            m = results[model]["test_metrics"]
            f1_thermo.append(m["f1_thermophile"])
            f1_nonthermo.append(m["f1_non_thermophile"])
            accuracies.append(m["accuracy"])

    x  = np.arange(len(models))
    w  = 0.22
    fig, ax = plt.subplots(figsize=(12, 5.5))

    def _split_yerr(model, key):
        if model == "THERMO":
            return None
        b  = results[model]["test_metrics_bootstrap"]
        v  = results[model]["test_metrics"][key]
        return [[v - b[key]["lower"]], [b[key]["upper"] - v]]

    thermo_errs    = [_split_yerr(m, "f1_thermophile")     for m in models]
    nonthermo_errs = [_split_yerr(m, "f1_non_thermophile") for m in models]
    acc_errs       = [_split_yerr(m, "accuracy")           for m in models]

    b1 = ax.bar(x - w, f1_thermo,    w, label="F1 — Thermophile",     color=CLASS_COLORS["thermophile"], edgecolor="white")
    b2 = ax.bar(x,     f1_nonthermo, w, label="F1 — Non-Thermophile", color=NON_EXTREMOPHILE, edgecolor="white")
    b3 = ax.bar(x + w, accuracies,   w, label="Accuracy",              color=BASELINE, edgecolor="white")

    for i, (bar1, bar2, bar3, e1, e2, e3) in enumerate(
            zip(b1, b2, b3, thermo_errs, nonthermo_errs, acc_errs)):
        if e1:
            ax.errorbar(bar1.get_x() + bar1.get_width()/2, f1_thermo[i],
                        yerr=[[e1[0][0]], [e1[1][0]]], fmt="none",
                        capsize=3, elinewidth=1.0, ecolor="#555", capthick=1.0)
        if e2:
            ax.errorbar(bar2.get_x() + bar2.get_width()/2, f1_nonthermo[i],
                        yerr=[[e2[0][0]], [e2[1][0]]], fmt="none",
                        capsize=3, elinewidth=1.0, ecolor="#555", capthick=1.0)
        if e3:
            ax.errorbar(bar3.get_x() + bar3.get_width()/2, accuracies[i],
                        yerr=[[e3[0][0]], [e3[1][0]]], fmt="none",
                        capsize=3, elinewidth=1.0, ecolor="#555", capthick=1.0)

    for bar in list(b1) + list(b2) + list(b3):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.008,
                f"{bar.get_height():.3f}", ha="center", va="bottom", fontsize=8)

    ax.set_xticks(x)
    ax.set_xticklabels([LABELS[m] for m in models], fontsize=10, rotation=12, ha="right")
    ax.set_ylim(0.5, 1.15)
    ax.set_ylabel("Score")
    ax.set_title("Binary Classifier — F1 per Class & Accuracy", fontsize=13)
    ax.legend(frameon=False, loc="upper left")

    plt.tight_layout()
    out = OUT_DIR / "binary_f1_breakdown.png"
    plt.savefig(out, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {out}")


if __name__ == "__main__":
    print("Generating binary classifier plots ...")
    plot_metric_comparison()
    plot_f1_breakdown()
    print(f"\nDone. All plots saved to {OUT_DIR}")
