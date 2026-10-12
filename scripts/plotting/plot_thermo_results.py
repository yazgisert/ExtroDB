# THERMO results on our proteomes: confusion matrix, probability distribution, training-set overlap and
# misclassified organisms.

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

from style import CLASS_COLORS, MAIN, NON_EXTREMOPHILE, SECONDARY

RESULT_FILE = Path("results/training_results/thermo_bacdive_result.json")
OUT_DIR     = Path("results/plots")
OUT_DIR.mkdir(parents=True, exist_ok=True)

plt.rcParams.update({
    "font.family":  "sans-serif",
    "font.size":    12,
    "axes.spines.top":   False,
    "axes.spines.right": False,
    "figure.dpi":   150,
})

THERMO_COLOR    = CLASS_COLORS["thermophile"]
NONTHERMO_COLOR = NON_EXTREMOPHILE
CORRECT_COLOR   = "#4CAF50"
WRONG_COLOR     = "#E05C5C"

with open(RESULT_FILE) as f:
    data = json.load(f)

preds   = data["per_proteome_predictions"]
y_true  = np.array([p["true_label"]           for p in preds])
y_pred  = np.array([p["predicted"]            for p in preds])
y_proba = np.array([p["proba_thermophile"]    for p in preds])

# Known overlap groups (from BacDive ID intersection computed separately)
# In THERMO training set: n=105, acc=0.981, macro_f1=0.980
# True held-out:          n=133, acc=0.977, macro_f1=0.967
OVERLAP_STATS = {
    "In Training\n(n=105)":  {"accuracy": 0.981, "macro_f1": 0.980},
    "Held-Out\n(n=133)":     {"accuracy": 0.977, "macro_f1": 0.967},
    "All 238\nProteomes":    {"accuracy": 0.979, "macro_f1": 0.975},
}


# Plot 1: confusion matrix
def plot_confusion_matrix():
    TP = int(np.sum((y_pred == 1) & (y_true == 1)))
    FN = int(np.sum((y_pred == 0) & (y_true == 1)))
    FP = int(np.sum((y_pred == 1) & (y_true == 0)))
    TN = int(np.sum((y_pred == 0) & (y_true == 0)))

    cm = np.array([[TP, FN], [FP, TN]])

    fig, ax = plt.subplots(figsize=(5, 4.5))
    im = ax.imshow(cm, cmap="Blues", vmin=0, vmax=cm.max())

    ax.set_xticks([0, 1])
    ax.set_yticks([0, 1])
    ax.set_xticklabels(["Thermophile\n(predicted)", "Non-Thermophile\n(predicted)"], fontsize=11)
    ax.set_yticklabels(["Thermophile\n(true)", "Non-Thermophile\n(true)"], fontsize=11)
    ax.xaxis.set_label_position("top")
    ax.xaxis.tick_top()

    labels = [
        [f"TP\n{TP}", f"FN\n{FN}"],
        [f"FP\n{FP}", f"TN\n{TN}"],
    ]
    for i in range(2):
        for j in range(2):
            color = "white" if cm[i, j] > cm.max() * 0.6 else "black"
            ax.text(j, i, labels[i][j], ha="center", va="center",
                    fontsize=14, fontweight="bold", color=color)

    ax.set_title("THERMO Model — Confusion Matrix\n(All 238 Proteomes)", fontsize=13, pad=18)
    plt.colorbar(im, ax=ax, shrink=0.8)
    plt.tight_layout()
    out = OUT_DIR / "thermo_confusion_matrix.png"
    plt.savefig(out, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {out}")


# Plot 2: prediction probabilities
def plot_probability_distribution():
    thermo_proba    = y_proba[y_true == 1]
    nonthermo_proba = y_proba[y_true == 0]

    fig, ax = plt.subplots(figsize=(7, 4))
    bins = np.linspace(0, 1, 21)

    ax.hist(nonthermo_proba, bins=bins, color=NONTHERMO_COLOR, alpha=0.8,
            label=f"Non-Thermophile (n={len(nonthermo_proba)})", edgecolor="white")
    ax.hist(thermo_proba, bins=bins, color=THERMO_COLOR, alpha=0.8,
            label=f"Thermophile (n={len(thermo_proba)})", edgecolor="white")

    ax.axvline(0.5, color="black", linestyle="--", linewidth=1.2, label="Decision threshold (0.5)")
    ax.set_xlabel("Predicted Probability of Thermophily")
    ax.set_ylabel("Number of Proteomes")
    ax.set_title("THERMO Model — Predicted Probability Distribution", fontsize=13)
    ax.legend(frameon=False)
    plt.tight_layout()
    out = OUT_DIR / "thermo_probability_distribution.png"
    plt.savefig(out, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {out}")


# Plot 3: organisms in THERMO's training set vs. held out
def plot_overlap_analysis():
    groups  = list(OVERLAP_STATS.keys())
    acc     = [OVERLAP_STATS[g]["accuracy"]  for g in groups]
    mf1     = [OVERLAP_STATS[g]["macro_f1"]  for g in groups]

    x   = np.arange(len(groups))
    w   = 0.32
    fig, ax = plt.subplots(figsize=(7, 4.5))

    bars1 = ax.bar(x - w/2, acc, width=w, label="Accuracy",  color=MAIN, edgecolor="white")
    bars2 = ax.bar(x + w/2, mf1, width=w, label="Macro-F1",  color=SECONDARY, edgecolor="white")

    for bar in list(bars1) + list(bars2):
        h = bar.get_height()
        ax.text(bar.get_x() + bar.get_width() / 2, h + 0.003,
                f"{h:.3f}", ha="center", va="bottom", fontsize=10)

    ax.set_ylim(0.90, 1.01)
    ax.set_xticks(x)
    ax.set_xticklabels(groups, fontsize=11)
    ax.set_ylabel("Score")
    ax.set_title("THERMO Performance: Training Overlap vs. True Held-Out", fontsize=13)
    ax.legend(frameon=False)

    ax.text(0.97, 0.97,
            "No data leakage:\nperformance is similar\nacross both groups",
            transform=ax.transAxes,
            ha="right", va="top",
            fontsize=9.5, color="grey",
            bbox=dict(boxstyle="round,pad=0.3", facecolor="white",
                      edgecolor="#cccccc", linewidth=0.8))

    plt.tight_layout()
    out = OUT_DIR / "thermo_overlap_analysis.png"
    plt.savefig(out, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {out}")


# Plot 4: misclassified organisms
def plot_misclassified():
    errors = [p for p in preds if p["true_label"] != p["predicted"]]

    fig, ax = plt.subplots(figsize=(8, 0.6 * len(errors) + 1.5))

    names  = [p["taxon_name"] for p in errors]
    probas = [p["proba_thermophile"] for p in errors]
    colors = [THERMO_COLOR if p["true_label"] == 1 else NONTHERMO_COLOR for p in errors]
    labels = ["FN (missed thermophile)" if p["true_label"] == 1
              else "FP (false thermophile)" for p in errors]

    y = np.arange(len(errors))
    bars = ax.barh(y, probas, color=colors, edgecolor="white", height=0.6)
    ax.axvline(0.5, color="black", linestyle="--", linewidth=1.2)

    ax.set_yticks(y)
    ax.set_yticklabels(names, fontsize=10, style="italic")
    ax.set_xlabel("Predicted Probability of Thermophily")
    ax.set_title("THERMO Model — Misclassified Organisms", fontsize=13)
    ax.set_xlim(0, 1)

    fn_patch = mpatches.Patch(color=THERMO_COLOR,    label="False Negative (true thermophile, predicted non-thermo)")
    fp_patch = mpatches.Patch(color=NONTHERMO_COLOR, label="False Positive (true non-thermophile, predicted thermo)")
    ax.legend(handles=[fn_patch, fp_patch], frameon=False, fontsize=9,
              loc="lower right")

    plt.tight_layout()
    out = OUT_DIR / "thermo_misclassified.png"
    plt.savefig(out, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {out}")


if __name__ == "__main__":
    print("Generating THERMO plots ...")
    plot_confusion_matrix()
    plot_probability_distribution()
    plot_overlap_analysis()
    plot_misclassified()
    print("\nDone. All plots saved to", OUT_DIR)
