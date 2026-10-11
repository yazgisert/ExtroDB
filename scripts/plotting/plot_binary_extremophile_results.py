# Plots for the binary extremophile classifier on the clustered split: all metrics and F1 per class, per embedder.

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

OUT_DIR = Path("results/plots")
OUT_DIR.mkdir(parents=True, exist_ok=True)

RESULT_FILES = {
    "One-Hot":  Path("results/training_results/binary_onehot_clustered_result.json"),
    "ProtT5":   Path("results/training_results/binary_prott5_clustered_result.json"),
    "ESM2":     Path("results/training_results/binary_esm2_650m_clustered_result.json"),
}

EMBEDDER_COLORS = {
    "One-Hot": "#A8A8A8",
    "ProtT5":  "#1B4F72",
    "ESM2":    "#E07840",
}

CLASS_COLORS = {
    "extremophile":     "#E05C5C",
    "non_extremophile": "#5C8AE0",
}

plt.rcParams.update({
    "font.family":       "sans-serif",
    "font.size":         12,
    "axes.spines.top":   False,
    "axes.spines.right": False,
    "figure.dpi":        150,
})

METRICS_KEYS   = ["accuracy", "balanced_accuracy", "macro_f1",
                   "f1_extremophile", "f1_non_extremophile", "mcc"]
METRICS_LABELS = ["Accuracy", "Balanced\nAccuracy", "Macro-F1",
                   "F1\n(Extremophile)", "F1\n(Non-Extrem.)", "MCC"]


def values_and_errors(result, keys):
    m     = result["test_metrics"]
    vals  = [m[k]["value"] for k in keys]
    lower = [m[k]["value"] - m[k]["ci_lower"] for k in keys]
    upper = [m[k]["ci_upper"] - m[k]["value"] for k in keys]
    return vals, lower, upper


def load_results():
    results = {}
    for name, path in RESULT_FILES.items():
        if path.exists():
            with open(path) as f:
                results[name] = json.load(f)
            print(f"  Loaded: {name} ({path.name})")
        else:
            print(f"  Skipping {name}, not found ({path.name})")
    return results


# Plot 1: all metrics, one embedder
def plot_single_metrics(result, embedder_name):
    vals, lower, upper = values_and_errors(result, METRICS_KEYS)

    fig, ax = plt.subplots(figsize=(11, 5))
    x    = np.arange(len(METRICS_KEYS))
    bars = ax.bar(x, vals, color=EMBEDDER_COLORS.get(embedder_name, "#1B4F72"),
                  edgecolor="white", width=0.55,
                  yerr=[lower, upper], capsize=4,
                  error_kw={"elinewidth": 1.2, "ecolor": "#666"})

    for bar, val in zip(bars, vals):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.015,
                f"{val:.3f}", ha="center", va="bottom", fontsize=10)

    ax.set_xticks(x)
    ax.set_xticklabels(METRICS_LABELS, fontsize=11)
    ax.set_ylim(0, 1.12)
    ax.set_ylabel("Score")
    ax.set_title(f"Binary Classifier — Extremophile vs. Non-Extremophile ({embedder_name})", fontsize=13)
    ax.axhline(1.0, color="grey", linewidth=0.5, linestyle="--", alpha=0.5)
    ax.text(0.99, 0.99,
            f"train: {result['n_train']:,} | val: {result['n_val']:,} | test: {result['n_test']:,}",
            transform=ax.transAxes, ha="right", va="top", fontsize=9, color="#888")

    plt.tight_layout()
    slug = embedder_name.lower().replace("-", "").replace(" ", "_")
    fname = f"binary_extremophile_metrics_{slug}_clustered.png"
    out = OUT_DIR / fname
    plt.savefig(out, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {out}")


# Plot 2: F1 per class and accuracy, one embedder
def plot_single_f1(result, embedder_name):
    keys   = ["f1_extremophile", "f1_non_extremophile", "accuracy"]
    labels = ["Extremophile", "Non-Extremophile", "Accuracy"]
    vals, lo, hi = values_and_errors(result, keys)
    colors = [CLASS_COLORS["extremophile"], CLASS_COLORS["non_extremophile"], "#888888"]
    fig, ax = plt.subplots(figsize=(8, 5))
    bars = ax.bar(labels, vals,
                  color=colors, edgecolor="white", width=0.4,
                  yerr=[lo, hi], capsize=4,
                  error_kw={"elinewidth": 1.2, "ecolor": "#333", "capthick": 1.2})
    for bar, val in zip(bars, vals):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.015,
                f"{val:.3f}", ha="center", va="bottom", fontsize=12, fontweight="bold")
    ax.set_ylim(0, 1.12)
    ax.set_ylabel("Score")
    ax.set_title(f"Binary Classifier — F1 per Class & Accuracy ({embedder_name})", fontsize=13)
    ax.axhline(1.0, color="grey", linewidth=0.5, linestyle="--", alpha=0.5)
    plt.tight_layout()
    slug  = embedder_name.lower().replace("-", "").replace(" ", "_")
    fname = f"binary_extremophile_f1_{slug}_clustered.png"
    out = OUT_DIR / fname
    plt.savefig(out, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {out}")


if __name__ == "__main__":
    print("Loading results ...")
    results = load_results()

    if not results:
        print("No result files found. Run training scripts first.")
    else:
        print("\nGenerating plots ...")
        for embedder, result in results.items():
            plot_single_metrics(result, embedder)
            plot_single_f1(result, embedder)
        print(f"\nDone. Plots saved to {OUT_DIR}")
