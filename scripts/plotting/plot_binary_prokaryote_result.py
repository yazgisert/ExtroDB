# Test-set metrics of the binary ProtT5 classifier on the clustered prokaryote split (n_test = 5,214).

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from style import MAIN, NON_EXTREMOPHILE, SECONDARY

OUT_DIR = Path("results/plots")
RESULT  = Path("results/training_results/binary_prott5_clustered_result.json")
OUT_DIR.mkdir(parents=True, exist_ok=True)

plt.rcParams.update({
    "font.family":       "sans-serif",
    "font.size":         12,
    "axes.spines.top":   False,
    "axes.spines.right": False,
    "figure.dpi":        150,
})


def main():
    with open(RESULT) as f:
        d = json.load(f)

    metrics = d["test_metrics"]

    items = [
        ("Accuracy",             "accuracy",           SECONDARY),
        ("Balanced Accuracy",    "balanced_accuracy",  SECONDARY),
        ("Macro F1",             "macro_f1",           SECONDARY),
        ("F1 — Extremophile",    "f1_extremophile",    MAIN),
        ("F1 — Non-extremophile","f1_non_extremophile",NON_EXTREMOPHILE),
        ("MCC",                  "mcc",                MAIN),
    ]

    labels = [it[0] for it in items]
    vals   = [metrics[it[1]]["value"] for it in items]
    colors = [it[2] for it in items]

    yerr_lo = [vals[i] - metrics[it[1]]["ci_lower"] for i, it in enumerate(items)]
    yerr_hi = [metrics[it[1]]["ci_upper"] - vals[i] for i, it in enumerate(items)]

    fig, ax = plt.subplots(figsize=(9, 5))
    x = np.arange(len(labels))

    bars = ax.bar(x, vals, color=colors, edgecolor="white", width=0.55,
                  yerr=[yerr_lo, yerr_hi], capsize=4,
                  error_kw={"ecolor": "#888", "elinewidth": 1.2})

    for bar, val in zip(bars, vals):
        ax.text(bar.get_x() + bar.get_width() / 2,
                bar.get_height() + 0.008,
                f"{val:.3f}",
                ha="center", va="bottom", fontsize=10, fontweight="bold")

    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=10, rotation=15, ha="right")
    ax.set_ylim(0, 1.05)
    ax.set_ylabel("Score")
    ax.axhline(1.0, color="#ccc", linewidth=0.8, linestyle="--")

    fig.text(0.99, -0.02,
             f"n_train = {d['n_train']:,} (50/50)  ·  n_test = {d['n_test']:,}  ·  non-extremophiles clustered at 30% (13,031 representatives)"
             f"     Error bars: 95% CI (bootstrap, 1000 iterations)",
             ha="right", va="top",
             fontsize=8.5, color="#888")

    ax.set_title(
        "Binary Classifier — Prokaryote Test Set\n"
        "Extremophile vs. Non-Extremophile (ProtT5 + FNN, clustered split)",
        fontsize=12
    )

    plt.tight_layout()
    out = OUT_DIR / "binary_prokaryote_result_clustered.png"
    plt.savefig(out, bbox_inches="tight")
    plt.close()
    print(f"Saved: {out}")


if __name__ == "__main__":
    main()
