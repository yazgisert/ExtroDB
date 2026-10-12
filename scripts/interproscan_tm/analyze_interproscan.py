# Pfam domain enrichment per class from the InterProScan results, for prokaryotes and algae, with heatmap and
# comparison plots.

import argparse
import csv
import json
import math
import sys
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "plotting"))
from style import CLASS_COLORS as CLASS_COLORS_STYLE

PROK_TSV  = Path("results/interproscan/prokaryote_results.tsv")
EUK_TSV   = Path("results/interproscan/eukaryote_results.tsv")
OUT_DIR   = Path("results/plots")
JSON_OUT  = Path("results/interproscan/domain_analysis.json")
LEAKED_IDS = Path("results/eukaryote_leakage/flagged_protein_ids.txt")
OUT_DIR.mkdir(parents=True, exist_ok=True)

CLASS_NAMES = ["acidophile", "alkaliphile", "halophile", "psychrophile", "thermophile"]

CLASS_COLORS = {
    "thermophile":  CLASS_COLORS_STYLE["thermophile"],
    "halophile":    CLASS_COLORS_STYLE["halophile"],
    "acidophile":   CLASS_COLORS_STYLE["acidophile"],
    "psychrophile": CLASS_COLORS_STYLE["psychrophile"],
    "alkaliphile":  CLASS_COLORS_STYLE["alkaliphile"],
}

plt.rcParams.update({
    "font.family":       "sans-serif",
    "font.size":         11,
    "axes.spines.top":   False,
    "axes.spines.right": False,
    "figure.dpi":        150,
})


def parse_pfam(tsv_path):
    pid2class   = {}
    pid2domains = defaultdict(set)
    domain_desc = {}

    with open(tsv_path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            if row["Analysis"] != "Pfam":
                continue
            raw   = row["Sequence_ID"]
            parts = raw.split("||")
            pid   = parts[0].split()[0]
            cls   = parts[-1].strip() if len(parts) > 1 else "unknown"

            pfam  = row["Signature_Accession"].strip()
            desc  = row["Signature_Description"].strip()

            pid2class[pid]   = cls
            pid2domains[pid].add(pfam)
            domain_desc[pfam] = desc

    return pid2class, pid2domains, domain_desc


def class_domain_matrix(pid2class, pid2domains, class_names):
    # Count total proteins per class (from all proteins seen, not just those with Pfam)
    class_protein_count = defaultdict(int)
    class_domain_counts = defaultdict(lambda: defaultdict(int))

    for pid, cls in pid2class.items():
        if cls not in class_names:
            continue
        class_protein_count[cls] += 1
        for pfam in pid2domains.get(pid, set()):
            class_domain_counts[cls][pfam] += 1

    return dict(class_protein_count), {k: dict(v) for k, v in class_domain_counts.items()}


def enrichment_score(class_freq, bg_freq, pseudocount=0.01):
    # log2 fold enrichment over the pooled background, with pseudocount
    return math.log2((class_freq + pseudocount) / (bg_freq + pseudocount))


def top_enriched(class_domain_counts, class_protein_count, class_names, top_n=15):
    # Background: pooled across all classes
    bg_counts = defaultdict(int)
    bg_total  = sum(class_protein_count.get(c, 0) for c in class_names)
    for cls in class_names:
        for pfam, cnt in class_domain_counts.get(cls, {}).items():
            bg_counts[pfam] += cnt

    results = {}
    for cls in class_names:
        n     = class_protein_count.get(cls, 1)
        dcnts = class_domain_counts.get(cls, {})
        bg_freq_all = {pfam: cnt / bg_total for pfam, cnt in bg_counts.items()}

        scored = []
        for pfam, cnt in dcnts.items():
            freq  = cnt / n
            bg    = bg_freq_all.get(pfam, 0)
            score = enrichment_score(freq, bg)
            scored.append((pfam, freq, score, cnt))

        scored.sort(key=lambda t: t[2], reverse=True)
        results[cls] = scored[:top_n]

    return results


# Plot 1: enrichment heatmap
def plot_heatmap(class_domain_counts, class_protein_count, domain_desc,
                 class_names, top_n=20, out_path=None):
    PSEUDOCOUNT = 0.005   # fraction pseudocount (0.5% of proteins)

    all_domains = set()
    for cls in class_names:
        all_domains.update(class_domain_counts.get(cls, {}).keys())

    # Per-domain mean frequency across classes (background)
    def mean_freq(pfam):
        freqs = [
            class_domain_counts.get(c, {}).get(pfam, 0) /
            max(class_protein_count.get(c, 1), 1)
            for c in class_names
        ]
        return sum(freqs) / len(freqs)

    def class_enrich(pfam, cls):
        n    = max(class_protein_count.get(cls, 1), 1)
        freq = class_domain_counts.get(cls, {}).get(pfam, 0) / n
        bg   = mean_freq(pfam)
        return math.log2((freq + PSEUDOCOUNT) / (bg + PSEUDOCOUNT))

    def max_enrich(pfam):
        return max(class_enrich(pfam, c) for c in class_names)

    # Select domains with at least one strong enrichment signal
    candidates = [d for d in all_domains if max_enrich(d) > 0.3]
    ranked = sorted(candidates, key=max_enrich, reverse=True)[:top_n]

    mat = np.zeros((len(ranked), len(class_names)))
    for j, cls in enumerate(class_names):
        for i, pfam in enumerate(ranked):
            mat[i, j] = class_enrich(pfam, cls)

    # Also build raw freq matrix for cell annotations
    raw = np.zeros_like(mat)
    for j, cls in enumerate(class_names):
        n = max(class_protein_count.get(cls, 1), 1)
        for i, pfam in enumerate(ranked):
            raw[i, j] = 100.0 * class_domain_counts.get(cls, {}).get(pfam, 0) / n

    row_labels = [f"{pfam}  {domain_desc.get(pfam, '')[:38]}" for pfam in ranked]
    col_labels = [c.capitalize() for c in class_names]

    fig_h = max(7, len(ranked) * 0.44)
    fig, ax = plt.subplots(figsize=(9, fig_h))

    vmax = max(abs(mat.max()), abs(mat.min()), 1.0)
    im = ax.imshow(mat, aspect="auto", cmap="RdBu_r", vmin=-vmax, vmax=vmax)
    cbar = plt.colorbar(im, ax=ax, shrink=0.55)
    cbar.set_label("log₂ enrichment vs. class mean\n(red = enriched, blue = depleted)", fontsize=9)

    ax.set_xticks(range(len(class_names)))
    ax.set_xticklabels(col_labels, fontsize=11, rotation=20, ha="right")
    ax.set_yticks(range(len(ranked)))
    ax.set_yticklabels(row_labels, fontsize=8)
    ax.set_title(
        "Pfam Domain Enrichment per Extremophile Class\n"
        "(Prokaryote sample · log₂ fold vs. class mean · raw % annotated)",
        fontsize=12,
    )

    for i in range(len(ranked)):
        for j in range(len(class_names)):
            r = raw[i, j]
            if r > 0:
                col = "white" if abs(mat[i, j]) > vmax * 0.6 else "black"
                ax.text(j, i, f"{r:.1f}%", ha="center", va="center",
                        fontsize=7, color=col)

    plt.tight_layout()
    plt.savefig(out_path, bbox_inches="tight")
    plt.close()
    print(f"Saved: {out_path}")


# Plot 2: top domains, prokaryotes vs. algae
def plot_prok_vs_euk(prok_counts, prok_n, euk_counts, euk_n,
                     domain_desc, class_names, top_n=8, out_path=None):
    n_cls = len(class_names)
    fig, axes = plt.subplots(1, n_cls, figsize=(4 * n_cls, 5.5), sharey=False)
    fig.suptitle("Top Pfam Domains: Prokaryote vs. Eukaryote\n(% proteins with domain, per class)",
                 fontsize=13, y=1.01)

    for ax, cls in zip(axes, class_names):
        pn = max(prok_n.get(cls, 1), 1)
        en = max(euk_n.get(cls, 1), 1)

        # Rank by prokaryote frequency
        p_dcnts = prok_counts.get(cls, {})
        ranked  = sorted(p_dcnts, key=lambda d: p_dcnts[d], reverse=True)[:top_n]

        if not ranked:
            ax.set_title(cls.capitalize(), fontsize=11)
            ax.axis("off")
            continue

        p_freqs = [100.0 * p_dcnts.get(d, 0) / pn for d in ranked]
        e_freqs = [100.0 * euk_counts.get(cls, {}).get(d, 0) / en for d in ranked]
        labels  = [f"{d}\n{domain_desc.get(d,'')[:22]}" for d in ranked]

        y = np.arange(len(ranked))
        h = 0.35

        color = CLASS_COLORS.get(cls, "#888")
        ax.barh(y + h/2, p_freqs, h, label="Prokaryote", color=color, alpha=0.9)
        ax.barh(y - h/2, e_freqs, h, label="Eukaryote",  color=color, alpha=0.40,
                hatch="//", edgecolor=color)

        ax.set_yticks(y)
        ax.set_yticklabels(labels, fontsize=7.5)
        ax.set_xlabel("% proteins", fontsize=9)
        ax.set_title(cls.capitalize(), fontsize=11, fontweight="bold", color=color)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        if ax == axes[0]:
            ax.legend(fontsize=8, frameon=False, loc="lower right")

    plt.tight_layout()
    plt.savefig(out_path, bbox_inches="tight")
    plt.close()
    print(f"Saved: {out_path}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--exclude-leaked", action="store_true",
                        help="drop algae proteins with >= 30%% identity to prokaryote train+val")
    args = parser.parse_args()
    suffix = "_filtered" if args.exclude_leaked else ""

    print("Parsing prokaryote results ...")
    p_pid2class, p_pid2domains, domain_desc = parse_pfam(PROK_TSV)
    p_class_n, p_class_dcnts = class_domain_matrix(p_pid2class, p_pid2domains, CLASS_NAMES)

    print("Parsing eukaryote results ...")
    e_pid2class, e_pid2domains, e_desc = parse_pfam(EUK_TSV)
    if args.exclude_leaked:
        leaked = {line.strip() for line in open(LEAKED_IDS)}
        n_before = len(e_pid2class)
        e_pid2class = {pid: cls for pid, cls in e_pid2class.items() if pid not in leaked}
        print(f"  Removed {n_before - len(e_pid2class)} leaked proteins, {len(e_pid2class)} left")
    domain_desc.update(e_desc)
    e_class_n, e_class_dcnts = class_domain_matrix(e_pid2class, e_pid2domains, CLASS_NAMES)

    print("\nProkaryote proteins with ≥1 Pfam hit per class:")
    for cls in CLASS_NAMES:
        print(f"  {cls:15s}: {p_class_n.get(cls, 0):4d} proteins, "
              f"{len(p_class_dcnts.get(cls, {})):3d} unique domains")

    print("\nEukaryote proteins with ≥1 Pfam hit per class:")
    for cls in CLASS_NAMES:
        print(f"  {cls:15s}: {e_class_n.get(cls, 0):4d} proteins, "
              f"{len(e_class_dcnts.get(cls, {})):3d} unique domains")

    prok_enriched = top_enriched(p_class_dcnts, p_class_n, CLASS_NAMES)
    print("\nTop 5 enriched Pfam domains per prokaryote class:")
    for cls in CLASS_NAMES:
        print(f"  {cls}:")
        for pfam, freq, score, cnt in prok_enriched[cls][:5]:
            print(f"    {pfam} ({domain_desc.get(pfam,'?')[:45]}) "
                  f"freq={freq:.3f} enrich={score:.2f}")

    print("\nGenerating plots ...")
    plot_heatmap(
        p_class_dcnts, p_class_n, domain_desc, CLASS_NAMES,
        top_n=20,
        out_path=OUT_DIR / "interproscan_heatmap_prok.png",
    )
    plot_prok_vs_euk(
        p_class_dcnts, p_class_n,
        e_class_dcnts, e_class_n,
        domain_desc, CLASS_NAMES,
        top_n=8,
        out_path=OUT_DIR / f"interproscan_prok_vs_euk{suffix}.png",
    )

    summary = {
        "prokaryote": {
            cls: {
                "n_proteins": p_class_n.get(cls, 0),
                "n_unique_domains": len(p_class_dcnts.get(cls, {})),
                "top_enriched": [
                    {"pfam": pf, "description": domain_desc.get(pf, ""),
                     "frequency": round(fr, 4), "log2_enrichment": round(sc, 3),
                     "count": ct}
                    for pf, fr, sc, ct in prok_enriched.get(cls, [])[:10]
                ],
            }
            for cls in CLASS_NAMES
        },
        "eukaryote": {
            cls: {
                "n_proteins": e_class_n.get(cls, 0),
                "n_unique_domains": len(e_class_dcnts.get(cls, {})),
            }
            for cls in CLASS_NAMES
        },
    }
    json_out = JSON_OUT.with_name(f"domain_analysis{suffix}.json")
    with open(json_out, "w") as f:
        json.dump(summary, f, indent=2)
    print(f"\nSaved JSON to {json_out}")
    print("Done.")


if __name__ == "__main__":
    main()
