# Embedder comparison plots: binary test metrics and algae recall for ProtT5, ESM2 650M, ESM2 8M and one-hot.

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

RESULTS = Path("results/training_results")
OUT_DIR = Path("results/plots")
OUT_DIR.mkdir(parents=True, exist_ok=True)

EMBEDDERS = {
    "onehot":    ("One-hot (20)",        "#A8A8A8"),
    "esm2_8m":   ("ESM2 8M (320)",       "#F2B38A"),
    "esm2_650m": ("ESM2 650M (1280)",    "#E07840"),
    "prott5":    ("ProtT5 (1024)",       "#1B4F72"),
}
PROKARYOTE_METRICS = {
    "accuracy":            "Accuracy",
    "balanced_accuracy":   "Balanced\naccuracy",
    "specificity":         "Specificity",
    "recall":              "Recall",
    "f1_extremophile":     "F1\nextremophile",
    "f1_non_extremophile": "F1 non-\nextremophile",
    "mcc":                 "MCC",
}
ALGAE_GROUPS = ["overall", "acidophile", "halophile", "psychrophile", "thermophile"]

plt.rcParams.update({
    "font.family":       "sans-serif",
    "font.size":         12,
    "axes.spines.top":   False,
    "axes.spines.right": False,
    "figure.dpi":        150,
})


def grouped_bars(ax, groups, series):
    width = 0.8 / len(series)
    x = np.arange(len(groups))
    for i, (label, color, values, lower, upper) in enumerate(series):
        offset = (i - (len(series) - 1) / 2) * width
        values, lower, upper = map(np.array, (values, lower, upper))
        ax.bar(x + offset, values, width, color=color, edgecolor="white", label=label,
               yerr=[values - lower, upper - values], capsize=3,
               error_kw={"elinewidth": 1.1, "ecolor": "#333"})
        for xi, val, top in zip(x + offset, values, upper):
            ax.text(xi, top + 0.01, f"{val:.2f}", ha="center", va="bottom", fontsize=8.5,
                    fontweight="bold" if label.startswith("ProtT5") else "normal")
    ax.set_xticks(x)
    ax.set_xticklabels(groups, fontsize=11)
    ax.set_ylim(0, 1.1)
    ax.legend(frameon=False, fontsize=10, loc="upper center", bbox_to_anchor=(0.5, -0.15), ncol=len(series))


def plot_prokaryote():
    series = []
    for key, (label, color) in EMBEDDERS.items():
        path = RESULTS / f"binary_{key}_clustered_result.json"
        if not path.exists():
            continue
        with open(path) as f:
            m = json.load(f)["test_metrics"]
        series.append((label, color,
                       [m[k]["value"] for k in PROKARYOTE_METRICS],
                       [m[k]["ci_lower"] for k in PROKARYOTE_METRICS],
                       [m[k]["ci_upper"] for k in PROKARYOTE_METRICS]))
        n_test = json.load(open(path))["n_test"]

    fig, ax = plt.subplots(figsize=(12, 5))
    grouped_bars(ax, list(PROKARYOTE_METRICS.values()), series)
    ax.set_ylabel("Score on prokaryote test set")
    ax.set_title(f"Binary Extremophile Classifier — Embedder Comparison (clustered split, n_test = {n_test:,})",
                 fontsize=12)
    plt.tight_layout()
    out = OUT_DIR / "binary_embedder_comparison_clustered.png"
    plt.savefig(out, bbox_inches="tight")
    plt.close()
    print(f"Saved: {out}")


def plot_algae():
    with open(RESULTS / "eukaryote_embedder_clustered_result.json") as f:
        result = json.load(f)

    series = []
    for key, (label, color) in EMBEDDERS.items():
        if key not in result["results"]:
            continue
        r = result["results"][key]
        rows = [r["overall"]] + [r["by_class"][c] for c in ALGAE_GROUPS[1:]]
        series.append((label, color,
                       [row["recall"] for row in rows],
                       [row["ci_lower"] for row in rows],
                       [row["ci_upper"] for row in rows]))

    fig, ax = plt.subplots(figsize=(11, 5))
    grouped_bars(ax, [g.capitalize() for g in ALGAE_GROUPS], series)
    ax.axvline(0.5, color="#ccc", linewidth=1)
    ax.set_ylabel("Recall on algae (fraction flagged as extremophile)")
    ax.set_title(f"Eukaryote Transfer by Embedder — Leakage-Filtered Algae (n = {result['n_algae_proteins']:,})",
                 fontsize=12)
    plt.tight_layout()
    out = OUT_DIR / "eukaryote_embedder_comparison_clustered.png"
    plt.savefig(out, bbox_inches="tight")
    plt.close()
    print(f"Saved: {out}")


if __name__ == "__main__":
    plot_prokaryote()
    plot_algae()
