# Dataset statistics plots: class distribution, train/val/test breakdown per class, and sequence length per class.

import re
from pathlib import Path
from collections import defaultdict

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

CLASSES = ["acidophile", "alkaliphile", "halophile", "psychrophile", "thermophile"]
COLORS  = ["#E05C5C", "#E0A85C", "#5C8AE0", "#7EC8A0", "#9B59B6"]

SPLITS = {
    "split1_20pct": {
        "train": "dataset_splits/rep_split1_20pct_train.faa",
        "val":   "dataset_splits/rep_split1_20pct_val.faa",
        "test":  "dataset_splits/rep_split1_20pct_test.faa",
    },
    "split2_20pct": {
        "train": "dataset_splits/rep_split2_20pct_train.faa",
        "val":   "dataset_splits/rep_split2_20pct_val.faa",
        "test":  "dataset_splits/rep_split2_20pct_test.faa",
    },
    "split1_30pct": {
        "train": "dataset_splits/rep_split1_train.faa",
        "val":   "dataset_splits/rep_split1_val.faa",
        "test":  "dataset_splits/rep_split1_test.faa",
    },
    "split2_30pct": {
        "train": "dataset_splits/rep_split2_train.faa",
        "val":   "dataset_splits/rep_split2_val.faa",
        "test":  "dataset_splits/rep_split2_test.faa",
    },
}

FULL_FAA = "dataset_splits/split_input_prokaryote_singlelabel.faa"


def parse_faa(path):
    records = []
    current = None
    seq_len = 0
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if current:
                    records.append((*current, seq_len))
                parts = line[1:].split("||")
                org   = parts[0].strip()
                prot  = parts[1].split()[0].strip()
                label = parts[2].strip()
                current = (org, prot, label)
                seq_len = 0
            else:
                seq_len += len(line)
        if current:
            records.append((*current, seq_len))
    return records  # (org, prot, label, seq_len)


# Plot 1: class distribution
def plot_class_distribution():
    # Use split1_20pct (train+val+test) as representative of the full dataset
    all_records = []
    for subset, path in SPLITS["split1_20pct"].items():
        all_records.extend(parse_faa(path))

    label_counts = defaultdict(int)
    org_per_class = defaultdict(set)
    for org, prot, label, _ in all_records:
        label_counts[label] += 1
        org_per_class[label].add(org)

    classes   = CLASSES
    prot_counts = [label_counts[c] for c in classes]
    org_counts  = [len(org_per_class[c]) for c in classes]

    x  = np.arange(len(classes))
    w  = 0.35
    fig, ax = plt.subplots(figsize=(9, 5))
    b1 = ax.bar(x - w/2, prot_counts, w, label="Proteins",  color="#5C8AE0", edgecolor="white")
    b2 = ax.bar(x + w/2, org_counts,  w, label="Organisms", color="#E0A85C", edgecolor="white")

    for bar in list(b1) + list(b2):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 200,
                f"{int(bar.get_height()):,}", ha="center", va="bottom", fontsize=9)

    ax.set_xticks(x)
    ax.set_xticklabels([c.capitalize() for c in classes], fontsize=11)
    ax.set_ylabel("Count")
    ax.set_title("Dataset Class Distribution — Proteins & Organisms", fontsize=13)
    ax.legend(frameon=False)
    plt.tight_layout()
    out = OUT_DIR / "dataset_class_distribution.png"
    plt.savefig(out, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {out}")


# Plot 2: train/val/test breakdown per class
def plot_split_breakdown():
    fig, axes = plt.subplots(2, 2, figsize=(13, 9), sharey=False)
    axes = axes.flatten()

    for ax, (split_name, paths) in zip(axes, SPLITS.items()):
        split_data = {}
        for subset, path in paths.items():
            records = parse_faa(path)
            counts  = defaultdict(int)
            for _, _, label, _ in records:
                counts[label] += 1
            split_data[subset] = [counts[c] for c in CLASSES]

        x   = np.arange(len(CLASSES))
        w   = 0.25
        ax.bar(x - w, split_data["train"], w, label="Train", color="#5C8AE0", edgecolor="white")
        ax.bar(x,     split_data["val"],   w, label="Val",   color="#7EC8A0", edgecolor="white")
        ax.bar(x + w, split_data["test"],  w, label="Test",  color="#E05C5C", edgecolor="white")

        ax.set_xticks(x)
        ax.set_xticklabels([c.capitalize() for c in CLASSES], fontsize=9, rotation=15)
        ax.set_title(split_name, fontsize=11)
        ax.set_ylabel("Proteins")
        ax.legend(frameon=False, fontsize=9)

    plt.suptitle("Train / Val / Test Protein Counts per Class", fontsize=13, y=1.01)
    plt.tight_layout()
    out = OUT_DIR / "dataset_split_breakdown.png"
    plt.savefig(out, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {out}")


# Plot 3: sequence length per class
def plot_seq_length_distribution():
    all_records = []
    for subset, path in SPLITS["split1_20pct"].items():
        all_records.extend(parse_faa(path))
    records = all_records
    lengths_per_class = defaultdict(list)
    for _, _, label, seq_len in records:
        lengths_per_class[label].append(seq_len)

    fig, ax = plt.subplots(figsize=(9, 5))
    for cls, color in zip(CLASSES, COLORS):
        lengths = lengths_per_class[cls]
        ax.hist(lengths, bins=50, alpha=0.6, color=color,
                label=f"{cls.capitalize()} (n={len(lengths):,})", edgecolor="none")

    ax.set_xlabel("Sequence Length (aa)")
    ax.set_ylabel("Number of Proteins")
    ax.set_title("Sequence Length Distribution per Class", fontsize=13)
    ax.set_xlim(0, 3000)
    ax.legend(frameon=False, fontsize=9)
    plt.tight_layout()
    out = OUT_DIR / "dataset_seq_length_distribution.png"
    plt.savefig(out, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {out}")


if __name__ == "__main__":
    print("Generating dataset statistics plots ...")
    plot_class_distribution()
    plot_split_breakdown()
    plot_seq_length_distribution()
    print(f"\nDone. All plots saved to {OUT_DIR}")
