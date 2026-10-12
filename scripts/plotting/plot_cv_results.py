# Cross-validation plots: accuracy and macro-F1 per fold, per-class F1, and CV mean vs. the original test result.

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from style import MAIN, SECONDARY

OUT_DIR = Path("results/plots")
OUT_DIR.mkdir(parents=True, exist_ok=True)

plt.rcParams.update({
    "font.family":       "sans-serif",
    "font.size":         12,
    "axes.spines.top":   False,
    "axes.spines.right": False,
    "figure.dpi":        150,
})

with open("results/training_results/split1_20pct_cv5_result.json") as f:
    cv20 = json.load(f)
with open("results/training_results/split1_30pct_cv5_result.json") as f:
    cv30 = json.load(f)

# Original single-split test results (for comparison)
ORIGINAL = {
    "20pct": {"accuracy": 0.7160, "macro-f1_score": 0.7147, "matthews-corr-coeff": 0.6462},
    "30pct": {"accuracy": 0.6730, "macro-f1_score": 0.6580, "matthews-corr-coeff": 0.6050},
}

CLASSES     = cv20["classes"]
COLOR_20    = SECONDARY
COLOR_30    = MAIN
COLOR_ORIG  = "#888888"


# Plot 1: accuracy and macro-F1 per fold
def plot_fold_variance():
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5), sharey=False)

    for ax, cv, color, label, orig_key in zip(
        axes,
        [cv20, cv30],
        [COLOR_20, COLOR_30],
        ["20% test split", "30% test split"],
        ["20pct", "30pct"],
    ):
        folds = list(range(1, 6))
        acc_folds = cv["cv_metrics_aggregated"]["accuracy"]["per_fold"]
        f1_folds  = cv["cv_metrics_aggregated"]["macro-f1_score"]["per_fold"]
        acc_mean  = cv["cv_metrics_aggregated"]["accuracy"]["mean"]
        f1_mean   = cv["cv_metrics_aggregated"]["macro-f1_score"]["mean"]
        acc_std   = cv["cv_metrics_aggregated"]["accuracy"]["std"]
        f1_std    = cv["cv_metrics_aggregated"]["macro-f1_score"]["std"]

        ax.plot(folds, acc_folds, "o-", color=color, linewidth=2,
                markersize=7, label="Accuracy (per fold)")
        ax.plot(folds, f1_folds, "s--", color=color, linewidth=2,
                markersize=7, alpha=0.7, label="Macro-F1 (per fold)")

        ax.axhline(acc_mean, color=color, linewidth=1.2, linestyle=":",
                   label=f"Accuracy mean ({acc_mean:.3f}±{acc_std:.3f})")
        ax.axhline(f1_mean,  color=color, linewidth=1.2, linestyle="-.",
                   label=f"Macro-F1 mean ({f1_mean:.3f}±{f1_std:.3f})")

        ax.axhline(ORIGINAL[orig_key]["accuracy"], color=COLOR_ORIG,
                   linewidth=1.5, linestyle="--",
                   label=f"Original test accuracy ({ORIGINAL[orig_key]['accuracy']:.3f})")

        ax.set_xticks(folds)
        ax.set_xticklabels([f"Fold {i}" for i in folds])
        ax.set_ylim(0.45, 0.85)
        ax.set_ylabel("Score")
        ax.set_title(f"5-Fold CV — {label}", fontsize=13)
        ax.legend(frameon=False, fontsize=9, loc="lower right")

    plt.tight_layout()
    out = OUT_DIR / "cv_fold_variance.png"
    plt.savefig(out, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {out}")


# Plot 2: per-class F1, mean ± std
def plot_per_class_f1():
    fig, ax = plt.subplots(figsize=(9, 5))

    x   = np.arange(len(CLASSES))
    w   = 0.3

    means_20 = [cv20["cv_per_class_aggregated"][c]["f1_score"]["mean"] for c in CLASSES]
    stds_20  = [cv20["cv_per_class_aggregated"][c]["f1_score"]["std"]  for c in CLASSES]
    means_30 = [cv30["cv_per_class_aggregated"][c]["f1_score"]["mean"] for c in CLASSES]
    stds_30  = [cv30["cv_per_class_aggregated"][c]["f1_score"]["std"]  for c in CLASSES]

    bars1 = ax.bar(x - w/2, means_20, w, yerr=stds_20,
                   label="20% test split", color=COLOR_20,
                   capsize=5, error_kw={"linewidth": 1.5}, edgecolor="white")
    bars2 = ax.bar(x + w/2, means_30, w, yerr=stds_30,
                   label="30% test split", color=COLOR_30,
                   capsize=5, error_kw={"linewidth": 1.5}, edgecolor="white")

    for bar, val in zip(list(bars1) + list(bars2),
                        means_20 + means_30):
        ax.text(bar.get_x() + bar.get_width() / 2,
                bar.get_height() + 0.01,
                f"{val:.2f}", ha="center", va="bottom", fontsize=9)

    ax.set_xticks(x)
    ax.set_xticklabels([c.capitalize() for c in CLASSES], fontsize=11)
    ax.set_ylim(0, 1.05)
    ax.set_ylabel("Mean F1-Score (± std across folds)")
    ax.set_title("5-Fold CV — Per-Class F1 by Split", fontsize=13)
    ax.legend(frameon=False)
    ax.axhline(0.2, color="grey", linewidth=0.8, linestyle="--", alpha=0.5,
               label="Random chance")

    plt.tight_layout()
    out = OUT_DIR / "cv_per_class_f1.png"
    plt.savefig(out, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {out}")


# Plot 3: CV mean vs. original test result
def plot_cv_vs_original():
    metrics_to_show = ["accuracy", "macro-f1_score", "matthews-corr-coeff"]
    labels          = ["Accuracy", "Macro-F1", "MCC"]

    x = np.arange(len(metrics_to_show))
    w = 0.22

    cv20_means = [cv20["cv_metrics_aggregated"][m]["mean"] for m in metrics_to_show]
    cv20_stds  = [cv20["cv_metrics_aggregated"][m]["std"]  for m in metrics_to_show]
    cv30_means = [cv30["cv_metrics_aggregated"][m]["mean"] for m in metrics_to_show]
    cv30_stds  = [cv30["cv_metrics_aggregated"][m]["std"]  for m in metrics_to_show]
    orig20     = [ORIGINAL["20pct"][m] for m in metrics_to_show]
    orig30     = [ORIGINAL["30pct"][m] for m in metrics_to_show]

    fig, ax = plt.subplots(figsize=(9, 5))

    ax.bar(x - 1.5*w, orig20,    w, label="Original test (20%)", color=COLOR_20,   alpha=0.5, edgecolor="white")
    ax.bar(x - 0.5*w, cv20_means, w, label="CV mean (20%)",       color=COLOR_20,
           yerr=cv20_stds, capsize=5, error_kw={"linewidth": 1.5}, edgecolor="white")
    ax.bar(x + 0.5*w, orig30,    w, label="Original test (30%)", color=COLOR_30,   alpha=0.5, edgecolor="white")
    ax.bar(x + 1.5*w, cv30_means, w, label="CV mean (30%)",       color=COLOR_30,
           yerr=cv30_stds, capsize=5, error_kw={"linewidth": 1.5}, edgecolor="white")

    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=12)
    ax.set_ylim(0, 0.85)
    ax.set_ylabel("Score")
    ax.set_title("Original Test Result vs. 5-Fold CV Mean", fontsize=13)
    ax.legend(frameon=False, fontsize=10)

    plt.tight_layout()
    out = OUT_DIR / "cv_vs_original.png"
    plt.savefig(out, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {out}")


if __name__ == "__main__":
    print("Generating CV plots ...")
    plot_fold_variance()
    plot_per_class_f1()
    plot_cv_vs_original()
    print(f"\nDone. All plots saved to {OUT_DIR}")
