# Algae overview plots: 5-class prokaryote vs. algae (overall and per-class F1) and binary algae recall per class and
# per species.

import json
import warnings
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from sklearn.metrics import accuracy_score, balanced_accuracy_score, f1_score, matthews_corrcoef

from style import CLASS_COLORS, MAIN, SECONDARY

RESULTS = Path("results/training_results")
OUT_DIR = Path("results/plots")

PROKARYOTE_JSON = RESULTS / "split1_30pct_result.json"
FIVECLASS_JSON  = RESULTS / "eukaryote_fiveclass_filtered_result.json"
BINARY_JSON     = RESULTS / "binary_eukaryote_eval_clustered_result.json"

SEED = 42
N_BOOTSTRAP = 1000
CLASS_COLORS = {
    "halophile":    CLASS_COLORS["halophile"],
    "psychrophile": CLASS_COLORS["psychrophile"],
    "acidophile":   CLASS_COLORS["acidophile"],
    "thermophile":  CLASS_COLORS["thermophile"],
}
PROKARYOTE_COLOR = MAIN
ALGAE_COLOR      = SECONDARY

plt.rcParams.update({
    "font.family":       "sans-serif",
    "font.size":         12,
    "axes.spines.top":   False,
    "axes.spines.right": False,
    "figure.dpi":        150,
})


def save(fig, name):
    fig.tight_layout()
    fig.savefig(OUT_DIR / name, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {OUT_DIR / name}")


def labels_from_confusion(cm):
    y_true, y_pred = [], []
    for i, row in enumerate(cm):
        for j, count in enumerate(row):
            y_true += [i] * count
            y_pred += [j] * count
    return np.array(y_true), np.array(y_pred)


def algae_overall(cm, present):
    # The confusion matrix holds every (true, predicted) pair, so it is enough to bootstrap the metrics
    y_true, y_pred = labels_from_confusion(cm)
    def metrics(t, p):
        return {
            "accuracy":          accuracy_score(t, p),
            "balanced-accuracy": balanced_accuracy_score(t, p),
            "macro-f1_score":    f1_score(t, p, labels=present, average="macro", zero_division=0),
            "matthews-corr-coeff": matthews_corrcoef(t, p),
        }
    point = metrics(y_true, y_pred)
    rng = np.random.default_rng(SEED)
    boot = defaultdict(list)
    for _ in range(N_BOOTSTRAP):
        idx = rng.choice(len(y_true), len(y_true), replace=True)
        for k, v in metrics(y_true[idx], y_pred[idx]).items():
            boot[k].append(v)
    return {k: (v, np.percentile(boot[k], 2.5), np.percentile(boot[k], 97.5)) for k, v in point.items()}


def wilson(p, n, z=1.96):
    centre = (p + z * z / (2 * n)) / (1 + z * z / n)
    half = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / (1 + z * z / n)
    return centre - half, centre + half


def plot_overall(prok, algae):
    keys = ["accuracy", "balanced-accuracy", "macro-f1_score", "matthews-corr-coeff"]
    names = ["Accuracy", "Balanced\nAccuracy", "Macro-F1", "MCC"]
    p_vals = [prok["test_metrics"][k] for k in keys]
    p_err = [[p_vals[i] - prok["test_metrics_bootstrap"][k]["lower"] for i, k in enumerate(keys)],
             [prok["test_metrics_bootstrap"][k]["upper"] - p_vals[i] for i, k in enumerate(keys)]]
    a_vals = [algae[k][0] for k in keys]
    a_err = [[algae[k][0] - algae[k][1] for k in keys], [algae[k][2] - algae[k][0] for k in keys]]

    x, w = np.arange(len(keys)), 0.36
    fig, ax = plt.subplots(figsize=(9, 5))
    ax.bar(x - w / 2, p_vals, w, color=PROKARYOTE_COLOR, yerr=p_err, capsize=4, label="Prokaryote test (n = 28,576)")
    ax.bar(x + w / 2, a_vals, w, color=ALGAE_COLOR, yerr=a_err, capsize=4, label="Leakage-filtered algae (n = 133,689)")
    for xi, v in zip(x - w / 2, p_vals):
        ax.text(xi, v + 0.025, f"{v:.3f}", ha="center", fontsize=10)
    for xi, v in zip(x + w / 2, a_vals):
        ax.text(xi, v + 0.025, f"{v:.3f}", ha="center", fontsize=10)
    ax.set_xticks(x)
    ax.set_xticklabels(names)
    ax.set_ylim(0, 1.05)
    ax.set_ylabel("Score")
    ax.legend(frameon=False, fontsize=10)
    ax.set_title("5-Class Model — Prokaryote Test vs. Leakage-Filtered Algae (ProtT5)", fontsize=12)
    save(fig, "fiveclass_prokaryote_vs_algae_overall_filtered.png")


def plot_per_class_f1(prok, algae_f1):
    classes = ["acidophile", "halophile", "psychrophile", "thermophile"]
    p_vals = [prok["per_class_metrics"][c]["f1_score"] for c in classes]
    p_err = [[p_vals[i] - prok["per_class_metrics_bootstrap"][c]["f1_score"]["lower"] for i, c in enumerate(classes)],
             [prok["per_class_metrics_bootstrap"][c]["f1_score"]["upper"] - p_vals[i] for i, c in enumerate(classes)]]
    a_vals = [algae_f1[c]["value"] for c in classes]
    a_err = [[algae_f1[c]["value"] - algae_f1[c]["ci_lower"] for c in classes],
             [algae_f1[c]["ci_upper"] - algae_f1[c]["value"] for c in classes]]
    colors = [CLASS_COLORS[c] for c in classes]

    x, w = np.arange(len(classes)), 0.36
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.bar(x - w / 2, p_vals, w, color=colors, yerr=p_err, capsize=4)
    ax.bar(x + w / 2, a_vals, w, color=colors, alpha=0.45, hatch="//", edgecolor="white", yerr=a_err, capsize=4)
    for xi, v in zip(x - w / 2, p_vals):
        ax.text(xi, v + 0.025, f"{v:.2f}", ha="center", fontsize=10)
    for xi, v in zip(x + w / 2, a_vals):
        ax.text(xi, v + 0.025, f"{v:.2f}", ha="center", fontsize=10)
    ax.bar(0, 0, color="#888", label="Prokaryote test")
    ax.bar(0, 0, color="#888", alpha=0.45, hatch="//", label="Leakage-filtered algae")
    ax.set_xticks(x)
    ax.set_xticklabels([c.capitalize() for c in classes])
    ax.set_ylim(0, 1.05)
    ax.set_ylabel("F1 score")
    ax.legend(frameon=False, fontsize=10)
    ax.set_title("5-Class Model — Per-Class F1, Prokaryote Test vs. Filtered Algae", fontsize=12)
    save(fig, "fiveclass_prokaryote_vs_algae_per_class_f1_filtered.png")


def plot_binary_by_class(binary):
    classes = ["halophile", "psychrophile", "acidophile", "thermophile"]
    by_class = binary["by_class"]
    vals = [by_class[c]["recall"] for c in classes]
    err = [[vals[i] - by_class[c]["recall_ci_lower"] for i, c in enumerate(classes)],
           [by_class[c]["recall_ci_upper"] - vals[i] for i, c in enumerate(classes)]]

    fig, ax = plt.subplots(figsize=(8, 5))
    x = np.arange(len(classes))
    ax.bar(x, vals, 0.55, color=[CLASS_COLORS[c] for c in classes], yerr=err, capsize=4)
    ax.scatter(x, [by_class[c]["mean_prob"] for c in classes], color="black", zorder=3, label="Mean predicted prob.")
    for xi, v, c in zip(x, vals, classes):
        ax.text(xi, max(v, by_class[c]["mean_prob"]) + 0.03, f"{v:.3f}", ha="center", fontsize=11, fontweight="bold")
    ax.set_xticks(x)
    ax.set_xticklabels([c.capitalize() for c in classes])
    ax.set_ylim(0, 1.1)
    ax.set_ylabel("Recall (fraction flagged as extremophile)")
    ax.legend(frameon=False, fontsize=10, loc="lower left")
    ax.set_title(f"Binary Classifier (clustered split) — Recall on Leakage-Filtered Algae (n = {binary['n_algae_proteins']:,})", fontsize=11)
    save(fig, "binary_algae_per_class_recall_clustered.png")


def plot_binary_by_species(binary):
    order = {"thermophile": 0, "psychrophile": 1, "halophile": 2, "acidophile": 3}
    species = sorted(binary["by_species"].values(), key=lambda s: (order[s["broad_class"]], -s["recall"]))
    y = np.arange(len(species))[::-1]
    recalls = np.array([s["recall"] for s in species])
    lows, highs = zip(*(wilson(s["recall"], s["n_proteins"]) for s in species))

    fig, ax = plt.subplots(figsize=(9, 6))
    ax.barh(y, recalls, 0.65, color=[CLASS_COLORS[s["broad_class"]] for s in species],
            xerr=[recalls - np.array(lows), np.array(highs) - recalls], capsize=3)
    ax.scatter([s["mean_prob"] for s in species], y, color="black", s=14, zorder=3, label="Mean predicted prob.")
    for yi, s in zip(y, species):
        ax.text(1.01, yi, f"{s['recall']:.3f}  (n={s['n_proteins']:,})", va="center", fontsize=8)
    ax.set_yticks(y)
    ax.set_yticklabels([s["species"] for s in species], fontsize=9)
    ax.set_xlim(0, 1.25)
    ax.set_xlabel("Recall (fraction flagged as extremophile) · error bars = 95% Wilson CI")
    handles = [plt.Rectangle((0, 0), 1, 1, color=c) for c in CLASS_COLORS.values()]
    ax.legend(handles + [ax.collections[-1]], [c.capitalize() for c in CLASS_COLORS] + ["Mean predicted prob."],
              frameon=False, fontsize=9, loc="upper center", bbox_to_anchor=(0.5, -0.12), ncol=5)
    ax.set_title("Binary Classifier (clustered split) — Recall by Species, Leakage-Filtered Algae", fontsize=12)
    save(fig, "binary_algae_per_species_recall_clustered.png")


if __name__ == "__main__":
    # Algae have no alkaliphiles, but the model still predicts them; sklearn warns about that
    warnings.filterwarnings("ignore", message="y_pred contains classes not in y_true")
    prok = json.load(open(PROKARYOTE_JSON))
    fiveclass = json.load(open(FIVECLASS_JSON))["eukaryote_eval"]
    binary = json.load(open(BINARY_JSON))

    cm = np.array(fiveclass["confusion_matrix"])
    present = [i for i, row in enumerate(cm) if row.sum() > 0]
    plot_overall(prok, algae_overall(cm, present))
    plot_per_class_f1(prok, fiveclass["per_class_f1"])
    plot_binary_by_class(binary)
    plot_binary_by_species(binary)
