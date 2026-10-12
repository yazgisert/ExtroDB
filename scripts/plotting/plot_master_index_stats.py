# Master index plots: organisms per class, proteome sizes per class, and OGT / pH / NaCl distributions.

import csv
from collections import defaultdict, Counter
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from style import CLASS_COLORS, MAIN, SECONDARY, apply

MASTER_INDEX = Path("searchable_resource/extremophile_master_index.tsv")
OUT_DIR      = Path("results/plots")
OUT_DIR.mkdir(parents=True, exist_ok=True)

apply()

with open(MASTER_INDEX) as f:
    rows = list(csv.DictReader(f, delimiter="\t"))

# Map broad class from phenotype_labels
CLASS_MAP = {
    "Thermophile":          "thermophile",
    "Hyperthermophile":     "thermophile",
    "Thermophilic Alga":    "thermophile",
    "Thermoacidophilic Alga": "thermophile",
    "High-Temperature Alga": "thermophile",
    "Psychrophile":         "psychrophile",
    "Psychrophilic Alga":   "psychrophile",
    "Halophile (Slight)":   "halophile",
    "Halophile (Moderate)": "halophile",
    "Halophile (Extreme)":  "halophile",
    "Halophilic Alga":      "halophile",
    "Saline-Adapted Alga":  "halophile",
    "Acidophile":           "acidophile",
    "Moderate Acidophile":  "acidophile",
    "Hyperacidophile":      "acidophile",
    "Acidotolerant":        "acidophile",
    "Acidophilic Alga":     "acidophile",
    "Acidotolerant":        "acidophile",
    "Thermoacidophilic Alga": "acidophile",
    "Alkaliphile":          "alkaliphile",
    "Extreme Alkaliphile":  "alkaliphile",
    "Alkalitolerant":       "alkaliphile",
}

CLASSES = ["thermophile", "psychrophile", "halophile", "acidophile", "alkaliphile"]
COLORS  = [CLASS_COLORS[c] for c in CLASSES]
COLOR_MAP = dict(zip(CLASSES, COLORS))


def get_broad_classes(row):
    labels = [l.strip() for l in row["phenotype_labels"].split("|")]
    broad  = set()
    for l in labels:
        if l in CLASS_MAP:
            broad.add(CLASS_MAP[l])
    return broad


# Plot 1: organisms per class
def plot_organism_counts():
    class_counts  = Counter()
    multi_label   = 0
    domain_counts = {"prokaryote": 0, "eukaryote": 0}

    for row in rows:
        broad = get_broad_classes(row)
        for c in broad:
            class_counts[c] += 1
        if len(broad) > 1:
            multi_label += 1
        domain_counts[row["domain_scope"]] += 1

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # Left: organism count per class
    ax = axes[0]
    counts = [class_counts[c] for c in CLASSES]
    bars = ax.bar(CLASSES, counts, color=COLORS, edgecolor="white", width=0.6)
    for bar, val in zip(bars, counts):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.5,
                str(val), ha="center", va="bottom", fontsize=11, fontweight="bold")
    ax.set_xticklabels([c.capitalize() for c in CLASSES], fontsize=11)
    ax.set_ylabel("Number of Organisms")
    ax.set_title(f"Organisms per Class\n(total: {len(rows)}, multi-label: {multi_label})", fontsize=12)

    # Right: prokaryote vs eukaryote
    ax2 = axes[1]
    ax2.pie(
        [domain_counts["prokaryote"], domain_counts["eukaryote"]],
        labels=[f"Prokaryote\n({domain_counts['prokaryote']})",
                f"Eukaryote\n({domain_counts['eukaryote']})"],
        colors=[MAIN, SECONDARY],
        autopct="%1.1f%%", startangle=90,
        textprops={"fontsize": 11},
        wedgeprops={"edgecolor": "white"}
    )
    ax2.set_title("Domain Scope", fontsize=12)

    plt.suptitle("Master Index — Organism Overview", fontsize=13, y=1.02)
    plt.tight_layout()
    out = OUT_DIR / "master_organism_counts.png"
    plt.savefig(out, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {out}")


# Plot 2: proteome sizes per class
def plot_protein_counts():
    protein_counts = defaultdict(list)
    for row in rows:
        broad = get_broad_classes(row)
        try:
            n = int(row["protein_count"])
        except (ValueError, KeyError):
            continue
        for c in broad:
            protein_counts[c].append(n)

    fig, ax = plt.subplots(figsize=(9, 5))
    data   = [protein_counts[c] for c in CLASSES]
    bp = ax.boxplot(data, patch_artist=True, notch=False,
                    medianprops={"color": "black", "linewidth": 2})
    for patch, color in zip(bp["boxes"], COLORS):
        patch.set_facecolor(color)
        patch.set_alpha(0.7)

    ax.set_xticks(range(1, len(CLASSES)+1))
    ax.set_xticklabels([c.capitalize() for c in CLASSES], fontsize=11)
    ax.set_ylabel("Protein Count per Proteome")
    ax.set_title("Proteome Size Distribution per Class", fontsize=13)
    ax.set_ylim(0, 12000)
    ax.text(0.99, 0.98, "Outliers above 12,000 not shown", transform=ax.transAxes,
            ha="right", va="top", fontsize=8, color="grey")

    # Annotate medians — placed above Q3 (top of box) to avoid overlap
    for i, d in enumerate(data, 1):
        med = np.median(d)
        q3  = np.percentile(d, 75)
        ax.text(i, min(q3, 10500) + 300, f"median\n{int(med):,}", ha="center", va="bottom",
                fontsize=8.5, color="black",
                bbox=dict(boxstyle="round,pad=0.2", facecolor="white", edgecolor="none", alpha=0.8))

    plt.tight_layout()
    out = OUT_DIR / "master_protein_counts.png"
    plt.savefig(out, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {out}")


# Plot 3: OGT, pH and NaCl distributions
def plot_environmental_params():
    ogt   = defaultdict(list)
    ph    = defaultdict(list)
    nacl  = defaultdict(list)

    for row in rows:
        broad = get_broad_classes(row)
        try:
            if row["ogt_celsius"]:
                v = float(row["ogt_celsius"])
                for c in broad: ogt[c].append(v)
        except ValueError: pass
        try:
            if row["ph_optimum"]:
                v = float(row["ph_optimum"])
                for c in broad: ph[c].append(v)
        except ValueError: pass
        try:
            if row["nacl_optimum_M"]:
                v = float(row["nacl_optimum_M"])
                for c in broad: nacl[c].append(v)
        except ValueError: pass

    fig, axes = plt.subplots(1, 3, figsize=(14, 5.5))

    handles_dict = {}  # collect one handle per class for the shared legend
    for ax, data_dict, xlabel, title in zip(
        axes,
        [ogt, ph, nacl],
        ["Temperature (°C)", "pH Optimum", "NaCl Optimum (M)"],
        ["Optimal Growth Temperature", "Optimal pH", "Optimal NaCl Concentration"]
    ):
        for cls, color in zip(CLASSES, COLORS):
            vals = data_dict[cls]
            if vals:
                patches = ax.hist(vals, bins=12, alpha=0.6, color=color,
                                  label=f"{cls.capitalize()} (n={len(vals)})",
                                  edgecolor="none")[2]
                if cls not in handles_dict:
                    import matplotlib.patches as mpatches
                    handles_dict[cls] = mpatches.Patch(color=color, alpha=0.6,
                                                       label=cls.capitalize())
        ax.set_xlabel(xlabel)
        ax.set_ylabel("Organisms")
        ax.set_title(title, fontsize=11)
    fig.legend(handles=list(handles_dict.values()),
               loc="lower center", ncol=5, frameon=False,
               fontsize=9, bbox_to_anchor=(0.5, -0.06))

    plt.suptitle("Environmental Parameter Distributions per Class", fontsize=13, y=1.02)
    plt.tight_layout()
    plt.subplots_adjust(bottom=0.18)
    out = OUT_DIR / "master_env_params.png"
    plt.savefig(out, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {out}")


if __name__ == "__main__":
    print("Generating master index statistics plots ...")
    plot_organism_counts()
    plot_protein_counts()
    plot_environmental_params()
    print(f"\nDone. All plots saved to {OUT_DIR}")
