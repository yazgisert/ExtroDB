# 5-class model on non-extremophile proteins: confidence compared with extremophiles, and which classes they get.

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
from scipy.stats import gaussian_kde

OUT_DIR = Path("results/plots")
RESULT  = Path("results/training_results/nonextremophile_eval_result.json")

OUT_DIR.mkdir(parents=True, exist_ok=True)

CLASS_COLORS = {
    "thermophile":  "#E05C5C",
    "psychrophile": "#5C8AE0",
    "halophile":    "#E0A040",
    "acidophile":   "#6ABF69",
    "alkaliphile":  "#9B59B6",
}
CLASS_ORDER = ["thermophile", "psychrophile", "halophile", "acidophile", "alkaliphile"]

plt.rcParams.update({
    "font.family":       "sans-serif",
    "font.size":         12,
    "axes.spines.top":   False,
    "axes.spines.right": False,
    "figure.dpi":        150,
})


def main():
    with open(RESULT) as f:
        data = json.load(f)

    ext_probs  = np.array(data["extremophile_max_prob"])
    nxt_probs  = np.array(data["nonextreme_max_prob"])
    pred_pct   = data["nonextreme_predicted_class_pct"]
    pred_cnt   = data["nonextreme_predicted_class_counts"]
    n_nxt      = data["n_test_nonextreme"]

    ext_mean = data["extremophile_mean_max_prob"]
    nxt_mean = data["nonextreme_mean_max_prob"]

    fig, (ax_kde, ax_bar) = plt.subplots(1, 2, figsize=(13, 5))

    # Left panel: confidence (max softmax probability)
    x_range = np.linspace(0.15, 1.0, 400)

    kde_ext = gaussian_kde(ext_probs, bw_method=0.15)
    kde_nxt = gaussian_kde(nxt_probs, bw_method=0.15)

    ax_kde.fill_between(x_range, kde_ext(x_range), alpha=0.35,
                         color="#E05C5C", label=f"Extremophile proteins  (mean={ext_mean:.2f})")
    ax_kde.plot(x_range, kde_ext(x_range), color="#E05C5C", linewidth=1.8)

    ax_kde.fill_between(x_range, kde_nxt(x_range), alpha=0.35,
                         color="#5C8AE0", label=f"Non-extremophile proteins  (mean={nxt_mean:.2f})")
    ax_kde.plot(x_range, kde_nxt(x_range), color="#5C8AE0", linewidth=1.8)

    ax_kde.axvline(ext_mean, color="#E05C5C", linewidth=1.2,
                   linestyle="--", alpha=0.8)
    ax_kde.axvline(nxt_mean, color="#5C8AE0", linewidth=1.2,
                   linestyle="--", alpha=0.8)

    # Gap annotation
    y_ann = ax_kde.get_ylim()[1] * 0.6 if ax_kde.get_ylim()[1] > 0 else 3.0
    ax_kde.annotate(
        "", xy=(ext_mean, y_ann * 0.55), xytext=(nxt_mean, y_ann * 0.55),
        arrowprops=dict(arrowstyle="<->", color="#555", lw=1.2),
    )
    ax_kde.text((ext_mean + nxt_mean) / 2, y_ann * 0.58,
                f"Δ = {ext_mean - nxt_mean:.2f}",
                ha="center", va="bottom", fontsize=10, color="#555")

    ax_kde.set_xlabel("Max softmax probability (model confidence)")
    ax_kde.set_ylabel("Density")
    ax_kde.set_title("Confidence Distribution", fontsize=12)
    ax_kde.legend(frameon=False, fontsize=9)

    # Right panel: predicted classes of the non-extremophiles
    classes = CLASS_ORDER
    pcts    = [float(pred_pct.get(c, 0)) for c in classes]
    cnts    = [int(pred_cnt.get(c, 0))   for c in classes]
    colors  = [CLASS_COLORS[c] for c in classes]

    order  = sorted(range(len(classes)), key=lambda i: pcts[i], reverse=True)
    classes = [classes[i] for i in order]
    pcts    = [pcts[i]   for i in order]
    cnts    = [cnts[i]   for i in order]
    colors  = [colors[i] for i in order]

    x = np.arange(len(classes))
    bars = ax_bar.bar(x, pcts, color=colors, edgecolor="white", width=0.6)

    for bar, pct, cnt in zip(bars, pcts, cnts):
        ax_bar.text(bar.get_x() + bar.get_width() / 2,
                    bar.get_height() + 0.5,
                    f"{pct:.1f}%\n({cnt:,})",
                    ha="center", va="bottom", fontsize=9)

    ax_bar.set_xticks(x)
    ax_bar.set_xticklabels([c.capitalize() for c in classes], fontsize=11)
    ax_bar.set_ylabel("% of non-extremophile proteins assigned")
    ax_bar.set_ylim(0, max(pcts) * 1.28)
    ax_bar.set_title(f"Predicted Class Distribution\n(n = {n_nxt:,} non-extremophile proteins)",
                     fontsize=12)
    ax_bar.axhline(100 / 5, color="#aaa", linewidth=1.0, linestyle=":",
                   label="Uniform baseline (20%)")
    ax_bar.legend(frameon=False, fontsize=9)

    fig.suptitle(
        "5-Class Model Evaluated on Non-Extremophile Proteins (Out-of-Distribution)\n"
        "Model trained only on extremophile proteins — must always predict one of 5 extremophile classes",
        fontsize=11, y=1.02
    )

    plt.tight_layout()
    out = OUT_DIR / "nonextremophile_eval.png"
    plt.savefig(out, bbox_inches="tight")
    plt.close()
    print(f"Saved: {out}")


if __name__ == "__main__":
    main()
