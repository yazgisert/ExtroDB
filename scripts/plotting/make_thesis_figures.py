# Figures and tables of the thesis Results chapter, in one place: no titles (the captions carry the statement), the
# same fonts and sizes everywhere, colours from style.py and 95% bootstrap CIs as error bars. Writes PDF and PNG
# figures and LaTeX tables to results/thesis_figures/. The figure plan with the source of every figure is
# results_figure_plan.md in the thesis folder.
# Run from the repository root: python scripts/plotting/make_thesis_figures.py [figure names ...]

import argparse
import csv
import json
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import style

RESULTS = Path("results/training_results")
OUT_DIR = Path("results/thesis_figures")
MASTER_INDEX = Path("searchable_resource/extremophile_master_index.tsv")
DOMAINS = Path("curation/source_tables/proteome_domain_ncbi.tsv")
ALGAE_LABELS = Path("embeddings/eukaryote_algae_labels.tsv")

CLASSES = ["thermophile", "psychrophile", "halophile", "acidophile", "alkaliphile"]
EMBEDDERS = [("onehot", "Composition"), ("esm2_8m", "ESM2 8M"), ("esm2_650m", "ESM2 650M"), ("prott5", "ProtT5")]
RESULT_SUFFIX = {"prott5": "", "esm2_650m": "_esm2", "esm2_8m": "_esm2_8m", "onehot": "_onehot"}

# Widths in inches for a single-column thesis layout (text width about 6.3 in)
FULL_WIDTH = 6.3


def setup():
    style.apply()
    plt.rcParams.update({
        "font.size": 9,
        "axes.labelsize": 9,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "legend.fontsize": 8,
        "legend.frameon": False,
        "savefig.dpi": 300,
        "savefig.bbox": "tight",
    })
    OUT_DIR.mkdir(parents=True, exist_ok=True)


def save(fig, name):
    for ext in ("pdf", "png"):
        fig.savefig(OUT_DIR / f"{name}.{ext}")
    plt.close(fig)
    print(f"  {OUT_DIR / name}.pdf / .png")


def label(cls):
    return cls.capitalize()


def error_bars(ax, x, values, lower, upper, color="#334155"):
    yerr = [np.array(values) - np.array(lower), np.array(upper) - np.array(values)]
    ax.errorbar(x, values, yerr=yerr, fmt="none", ecolor=color, elinewidth=0.8, capsize=2)


# ── R0: dataset ─────────────────────────────────────────────────────────────────────────────────────────────────────
def load_dataset_rows():
    domain = {r["proteome_id"]: r["domain"] for r in csv.DictReader(open(DOMAINS), delimiter="\t")}
    algae_class = {}
    with open(ALGAE_LABELS) as f:
        next(f)
        for line in f:
            _, proteome, cls = line.rstrip("\n").split("\t")
            algae_class[proteome] = cls
    rows = [r for r in csv.DictReader(open(MASTER_INDEX), delimiter="\t") if r["proteome_status"] == "downloaded"]
    for r in rows:
        r["classes"] = r["phenotype_focus"].split("|")
        if r["domain_scope"] == "eukaryote":
            r["group"] = "algae"
            r["classes"] = [algae_class.get(r["proteome_id"], r["classes"][0])]
        elif len(r["classes"]) > 1:
            r["group"] = "multi"
        else:
            r["group"] = domain.get(r["proteome_id"], "bacteria")
    return rows


def fig_dataset_organisms():
    # Organisms per class: single-label archaea and bacteria (the organisms of the splits), prokaryotes with more than
    # one label (counted in every class they belong to) and algae (test set only)
    rows = load_dataset_rows()
    groups = [("archaea", "Archaea (single label)", style.RAMP[3]),
              ("bacteria", "Bacteria (single label)", style.RAMP[1]),
              ("multi", "Prokaryotes with several labels", style.RAMP[0]),
              ("algae", "Algae (test set)", style.BASELINE)]
    counts = {g: Counter() for g, _, _ in groups}
    for r in rows:
        for cls in r["classes"]:
            counts[r["group"]][cls] += 1

    fig, ax = plt.subplots(figsize=(FULL_WIDTH, 2.6))
    x = np.arange(len(CLASSES))
    bottom = np.zeros(len(CLASSES))
    for g, name, color in groups:
        values = np.array([counts[g][c] for c in CLASSES])
        ax.bar(x, values, 0.6, bottom=bottom, color=color, label=name, edgecolor="white", linewidth=0.5)
        for xi, v, b in zip(x, values, bottom):
            if v >= 4:
                ax.text(xi, b + v / 2, str(v), ha="center", va="center", fontsize=7,
                        color="white" if g in ("archaea", "bacteria") else "#1e293b")
        bottom += values
    for xi, total in zip(x, bottom):
        ax.text(xi, total + 1, str(int(total)), ha="center", va="bottom", fontsize=8)
    ax.set_xticks(x, [label(c) for c in CLASSES])
    ax.set_ylabel("Organisms")
    ax.set_ylim(0, bottom.max() * 1.12)
    ax.legend(loc="upper left", bbox_to_anchor=(1.0, 1.0))
    save(fig, "fig_r1_dataset_organisms")
    n_single = sum(1 for r in rows if r["group"] in ("archaea", "bacteria"))
    n_multi = sum(1 for r in rows if r["group"] == "multi")
    n_algae = sum(1 for r in rows if r["group"] == "algae")
    print(f"    {len(rows)} organisms: {n_single} single-label prokaryotes, {n_multi} with several labels, {n_algae} algae")


def fig_dataset_optima():
    # Reported growth optima of the single-label prokaryotes, one panel per parameter, with the class thresholds
    rows = [r for r in load_dataset_rows() if r["group"] in ("archaea", "bacteria")]
    panels = [("ogt_celsius", "OGT (°C)", ["psychrophile", "thermophile"], [15, 45, 80]),
              ("nacl_optimum_M", "NaCl optimum (M)", ["halophile"], [0.2, 0.85, 3.4]),
              ("ph_optimum", "pH optimum", ["acidophile", "alkaliphile"], [3.0, 5.0, 7.5, 9.0, 10.0])]
    fig, axes = plt.subplots(1, 3, figsize=(FULL_WIDTH, 2.4), gridspec_kw={"width_ratios": [2, 1, 2]})
    rng = np.random.default_rng(0)
    for ax, (field, ylabel, classes, thresholds) in zip(axes, panels):
        for i, cls in enumerate(classes):
            for dom, marker in (("bacteria", "o"), ("archaea", "^")):
                vals = [float(r[field]) for r in rows if cls in r["classes"] and r["group"] == dom and r[field]]
                xs = i + rng.uniform(-0.18, 0.18, len(vals))
                ax.scatter(xs, vals, s=12, marker=marker, color=style.CLASS_COLORS[cls], alpha=0.8,
                           edgecolors="white", linewidths=0.3)
        for t in thresholds:
            ax.axhline(t, color="#94a3b8", linewidth=0.6, linestyle=":", zorder=0)
        ax.set_xticks(range(len(classes)), [label(c) for c in classes])
        ax.set_xlim(-0.6, len(classes) - 0.4)
        ax.set_ylabel(ylabel)
    handles = [plt.Line2D([], [], marker="o", linestyle="", color="#64748b", label="Bacteria"),
               plt.Line2D([], [], marker="^", linestyle="", color="#64748b", label="Archaea")]
    fig.tight_layout()
    fig.legend(handles=handles, loc="lower center", ncol=2, bbox_to_anchor=(0.5, 1.0))
    save(fig, "fig_r2_dataset_optima")
    missing = Counter(c for r in rows for c in r["classes"]
                      if not r[{"thermophile": "ogt_celsius", "psychrophile": "ogt_celsius", "halophile": "nacl_optimum_M",
                                 "acidophile": "ph_optimum", "alkaliphile": "ph_optimum"}[c]])
    if missing:
        print(f"    organisms without a value for their class parameter: {dict(missing)}")


def fig_proteome_sizes():
    # Proteome sizes of the single-label prokaryotes per class, with the class median
    rows = [r for r in load_dataset_rows() if r["group"] in ("archaea", "bacteria")]
    fig, ax = plt.subplots(figsize=(FULL_WIDTH, 2.4))
    rng = np.random.default_rng(0)
    for i, cls in enumerate(CLASSES):
        sizes = [int(r["protein_count"]) for r in rows if r["classes"] == [cls]]
        for dom, marker in (("bacteria", "o"), ("archaea", "^")):
            vals = [int(r["protein_count"]) for r in rows if r["classes"] == [cls] and r["group"] == dom]
            ax.scatter(i + rng.uniform(-0.18, 0.18, len(vals)), vals, s=12, marker=marker,
                       color=style.CLASS_COLORS[cls], alpha=0.8, edgecolors="white", linewidths=0.3)
        ax.hlines(np.median(sizes), i - 0.28, i + 0.28, color="#334155", linewidth=1.2)
    ax.set_xticks(range(len(CLASSES)), [label(c) for c in CLASSES])
    ax.set_xlim(-0.6, len(CLASSES) - 0.4)
    ax.set_ylabel("Proteins per proteome")
    handles = [plt.Line2D([], [], marker="o", linestyle="", color="#64748b", label="Bacteria"),
               plt.Line2D([], [], marker="^", linestyle="", color="#64748b", label="Archaea"),
               plt.Line2D([], [], color="#334155", label="Median")]
    ax.legend(handles=handles, loc="upper left", bbox_to_anchor=(1.0, 1.0))
    save(fig, "fig_data_proteome_sizes")
    sizes = [int(r["protein_count"]) for r in rows]
    print(f"    {len(rows)} proteomes, {sum(sizes):,} proteins, {min(sizes):,} to {max(sizes):,} (median {np.median(sizes):,.0f})")


# ── R1: classification ──────────────────────────────────────────────────────────────────────────────────────────────
def multiclass_result(split, embedder):
    return json.load(open(RESULTS / f"{split}{RESULT_SUFFIX[embedder]}_result.json"))


def fig_fiveclass_per_class():
    # Per-class F1 of the five-class models, 30% clustering, all four representations
    fig, ax = plt.subplots(figsize=(FULL_WIDTH, 2.6))
    x = np.arange(len(CLASSES))
    w = 0.2
    for i, (key, name) in enumerate(EMBEDDERS):
        r = multiclass_result("split1_30pct", key)
        vals = [r["per_class_metrics"][c]["f1_score"] for c in CLASSES]
        ci = [r["per_class_metrics_bootstrap"][c]["f1_score"] for c in CLASSES]
        xs = x + (i - 1.5) * w
        ax.bar(xs, vals, w, color=style.EMBEDDER_COLORS[key], label=name)
        error_bars(ax, xs, vals, [c["lower"] for c in ci], [c["upper"] for c in ci])
    ax.set_xticks(x, [label(c) for c in CLASSES])
    ax.set_ylabel("F1 score")
    ax.set_ylim(0, 1)
    ax.legend(ncol=4, loc="upper center", bbox_to_anchor=(0.5, 1.12))
    save(fig, "fig_r3_fiveclass_per_class_f1")


def fmt_ci(value, ci):
    return f"{value:.3f} {{\\scriptsize[{ci['lower']:.3f}, {ci['upper']:.3f}]}}"


def table_multiclass(pct="30pct", name="table_r1_multiclass"):
    # LaTeX table: five-class and twelve-class models, all representations, accuracy / macro-F1 / MCC with 95% CIs
    metrics = [("accuracy", "Accuracy"), ("macro-f1_score", "Macro-F1"), ("matthews-corr-coeff", "MCC")]
    lines = [r"\begin{tabular}{llccc}", r"\hline",
             r"\textbf{Task} & \textbf{Representation} & " + " & ".join(rf"\textbf{{{m}}}" for _, m in metrics) + r" \\",
             r"\hline"]
    for split, task in ((f"split1_{pct}", "Five-class"), (f"split2_{pct}", "Twelve-class")):
        for i, (key, rep_name) in enumerate(reversed(EMBEDDERS)):
            r = multiclass_result(split, key)
            cells = [fmt_ci(r["test_metrics"][m], r["test_metrics_bootstrap"][m]) for m, _ in metrics]
            lines.append(f"{task if i == 0 else ''} & {rep_name} & " + " & ".join(cells) + r" \\")
        lines.append(r"\hline")
    lines.append(r"\end{tabular}")
    path = OUT_DIR / f"{name}.tex"
    path.write_text("\n".join(lines) + "\n")
    print(f"  {path}")



# ── R2: robustness ──────────────────────────────────────────────────────────────────────────────────────────────────
LEAKAGE_SUMMARIES = {"initial": Path("results/test_vs_train_search/leakage_summary.txt"),
                     "clustered": Path("results/test_vs_train_search_clustered/leakage_summary.txt")}


def leaked_counts(path):
    # (leaked, total) per group, parsed from the summary written by analyze_test_leakage.py
    counts, group = {}, None
    for line in path.read_text().splitlines():
        if line.startswith(("extremophile:", "non_extremophile:")):
            group = line.split(":")[0]
            counts[group] = [None, int(line.split("/")[1].split("(")[0].replace(",", "").strip())]
        elif "leaked" in line and group:
            counts[group][0] = int(line.split(":")[1].split("(")[0].replace(",", "").strip())
    return counts


def table_leakage():
    # Share of leaked test proteins and binary MCC on the initial random split and on the clustered split
    leak = {k: leaked_counts(p) for k, p in LEAKAGE_SUMMARIES.items()}
    lines = [r"\begin{tabular}{lcc}", r"\hline",
             r" & \textbf{Initial random split} & \textbf{Clustered split} \\", r"\hline",
             r"\multicolumn{3}{l}{\textit{Leaked test proteins}} \\"]
    for group, name in (("non_extremophile", "Non-extremophiles"), ("extremophile", "Extremophiles")):
        cells = []
        for k in ("initial", "clustered"):
            n, total = leak[k][group]
            cells.append(f"{n:,} of {total:,} ({100 * n / total:.1f}\\%)".replace(",", "{,}"))
        lines.append(f"{name} & " + " & ".join(cells) + r" \\")
    lines += [r"\hline", r"\multicolumn{3}{l}{\textit{MCC of the binary classifier}} \\"]
    for key, name in reversed(EMBEDDERS):
        a = json.load(open(RESULTS / f"binary_{key}_samesplit_result.json"))["test_metrics"]["mcc"]
        b = json.load(open(RESULTS / f"binary_{key}_clustered_result.json"))["test_metrics"]["mcc"]
        cell = lambda m: f"{m['value']:.3f} {{\\scriptsize[{m['ci_lower']:.3f}, {m['ci_upper']:.3f}]}}"
        lines.append(f"{name} & {cell(a)} & {cell(b)}" + r" \\")
    lines += [r"\hline", r"\end{tabular}"]
    path = OUT_DIR / "table_r2_leakage.tex"
    path.write_text("\n".join(lines) + "\n")
    print(f"  {path}")


def fig_homology():
    # Homology-based inference vs. ProtT5 on all test proteins and on those with and without a hit in training
    fig, axes = plt.subplots(1, 2, figsize=(FULL_WIDTH, 2.6), sharey=True)
    for ax, (file, title) in zip(axes, (("homology_fiveclass_result.json", "Five-class test set"),
                                         ("homology_binary_clustered_result.json", "Binary test set"))):
        r = json.load(open(RESULTS / file))
        sub = r["subsets"]
        groups = [f"All\n{r['n_test']:,}", f"With hit\n{sub['with_hit']['n']:,}",
                  f"Without hit\n{sub['without_hit']['n']:,}"]
        hom = [(r["homology"]["accuracy"]["value"], r["homology"]["accuracy"]["ci_lower"], r["homology"]["accuracy"]["ci_upper"]),
               (sub["with_hit"]["homology_accuracy"], *sub["with_hit"]["homology_accuracy_ci"]),
               None]
        pt5 = [(r["prott5"]["accuracy"]["value"], r["prott5"]["accuracy"]["ci_lower"], r["prott5"]["accuracy"]["ci_upper"]),
               (sub["with_hit"]["prott5_accuracy"], *sub["with_hit"]["prott5_accuracy_ci"]),
               (sub["without_hit"]["prott5_accuracy"], *sub["without_hit"]["prott5_accuracy_ci"])]
        x = np.arange(3)
        w = 0.36
        for offset, data, color, name in ((-w / 2, hom, style.BASELINE, "Homology-based inference"),
                                          (w / 2, pt5, style.MAIN, "ProtT5")):
            idx = [i for i, d in enumerate(data) if d]
            vals = [data[i][0] for i in idx]
            ax.bar(x[idx] + offset, vals, w, color=color, label=name)
            error_bars(ax, x[idx] + offset, vals, [data[i][1] for i in idx], [data[i][2] for i in idx])
        ax.text(x[2] - w / 2, 0.03, "no prediction", ha="center", va="bottom", rotation=90, fontsize=7, color="#1e293b")
        ax.set_xticks(x, groups, fontsize=7.5)
        ax.set_title(title, fontsize=9)
        ax.set_ylim(0, 1)
    axes[0].set_ylabel("Accuracy")
    handles, labels_ = axes[0].get_legend_handles_labels()
    fig.tight_layout(w_pad=2)
    fig.legend(handles, labels_, loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=2)
    save(fig, "fig_r4_homology_vs_prott5")


def fig_crossval():
    # Per-protein split vs. organism-grouped cross-validation (mean and standard deviation over five folds)
    fig, ax = plt.subplots(figsize=(FULL_WIDTH * 0.62, 2.5))
    metrics = [("accuracy", "Accuracy"), ("macro-f1_score", "Macro-F1")]
    labels, x, w = [], np.arange(4), 0.36
    split_vals, split_lo, split_hi, cv_mean, cv_sd = [], [], [], [], []
    for pct in ("30pct", "20pct"):
        r = multiclass_result(f"split1_{pct}", "prott5")
        cv = json.load(open(RESULTS / f"split1_{pct}_cv5_result.json"))["cv_metrics_aggregated"]
        for m, name in metrics:
            labels.append(f"{name}\n({pct[:2]}%)")
            split_vals.append(r["test_metrics"][m])
            split_lo.append(r["test_metrics_bootstrap"][m]["lower"])
            split_hi.append(r["test_metrics_bootstrap"][m]["upper"])
            cv_mean.append(cv[m]["mean"])
            cv_sd.append(cv[m]["std"])
    ax.bar(x - w / 2, split_vals, w, color=style.MAIN, label="Per-protein split (95% CI)")
    error_bars(ax, x - w / 2, split_vals, split_lo, split_hi)
    ax.bar(x + w / 2, cv_mean, w, color=style.SECONDARY, label="Unseen organisms, 5-fold CV (mean ± SD)")
    ax.errorbar(x + w / 2, cv_mean, yerr=cv_sd, fmt="none", ecolor="#334155", elinewidth=0.8, capsize=2)
    ax.axhline(0.2, color="#94a3b8", linewidth=0.7, linestyle=":")
    ax.set_xlim(-0.6, 4.1)
    ax.text(3.62, 0.2, "chance", fontsize=7, color="#1e293b", ha="left", va="center")
    ax.set_xticks(x, labels, fontsize=7.5)
    ax.set_ylabel("Score")
    ax.set_ylim(0, 1)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, 1.22), ncol=1)
    save(fig, "fig_r5_crossval")


def fig_archaea():
    # A: five-class ProtT5 accuracy per class for archaeal and bacterial test proteins.
    # B: halophile and thermophile subtypes by domain, ProtT5 vs. composition (single-organism groups left out)
    pt5 = json.load(open(RESULTS / "fiveclass_domain_bias_result.json"))
    comp = json.load(open(RESULTS / "fiveclass_domain_bias_onehot_result.json"))
    fig, (ax_a, ax_b) = plt.subplots(2, 1, figsize=(FULL_WIDTH, 4.8))

    x, w = np.arange(len(CLASSES)), 0.38
    for offset, dom, color, name in ((-w / 2, "archaea", style.RAMP[3], "Archaeal proteins"),
                                     (w / 2, "bacteria", style.RAMP[1], "Bacterial proteins")):
        for xi, cls in zip(x, CLASSES):
            d = pt5["per_class"][cls].get(dom)
            if d and d["n"] >= 100:
                ax_a.bar(xi + offset, d["accuracy"], w, color=color, label=name if cls == "halophile" else None)
                ax_a.text(xi + offset, d["accuracy"] + 0.02, f"{d['accuracy']:.2f}", ha="center", fontsize=6)
            elif not d or d["n"] == 0:
                ax_a.text(xi + offset, 0.03, "none", ha="center", rotation=90, fontsize=7, color="#1e293b")
    ax_a.set_xticks(x, [label(c) for c in CLASSES])
    ax_a.set_ylabel("Accuracy (= recall)")
    ax_a.set_ylim(0, 1.08)
    ax_a.legend(loc="upper left", bbox_to_anchor=(1.0, 1.0))
    ax_a.text(-0.1, 1.02, "A", transform=ax_a.transAxes, fontweight="bold", fontsize=10)

    groups = [("extreme_halophile|archaea", "Extreme\narchaea"), ("moderate_halophile|archaea", "Moderate\narchaea"),
              ("moderate_halophile|bacteria", "Moderate\nbacteria"), ("slight_halophile|bacteria", "Slight\nbacteria"),
              ("hyperthermophile|archaea", "Hyper-\nthermophile\narchaea"), ("thermophile|archaea", "Thermo-\nphile\narchaea"),
              ("thermophile|bacteria", "Thermo-\nphile\nbacteria")]
    xb = np.arange(len(groups))
    for offset, res, color, name in ((-w / 2, pt5, style.MAIN, "ProtT5"), (w / 2, comp, style.RAMP[0], "Composition")):
        vals = [res["accuracy_by_subtype_and_domain"][g]["accuracy"] for g, _ in groups]
        ax_b.bar(xb + offset, vals, w, color=color, label=name)
        for xi, v in zip(xb + offset, vals):
            ax_b.text(xi, v + 0.02, f"{v:.2f}", ha="center", fontsize=6)
    ax_b.axvline(3.5, color="#94a3b8", linewidth=0.7, linestyle="--")
    ax_b.text(1.5, 1.04, "Halophile subtypes", ha="center", fontsize=7.5, color="#475569")
    ax_b.text(5.0, 1.04, "Thermophile subtypes", ha="center", fontsize=7.5, color="#475569")
    ax_b.set_xticks(xb, [n for _, n in groups], fontsize=7.5)
    ax_b.set_ylabel("Accuracy (= recall)")
    ax_b.set_ylim(0, 1.1)
    ax_b.legend(loc="upper left", bbox_to_anchor=(1.0, 1.0))
    ax_b.text(-0.1, 1.02, "B", transform=ax_b.transAxes, fontweight="bold", fontsize=10)
    fig.tight_layout()
    save(fig, "fig_r6_archaea")


def pct_ci(v, digits=1, scale=100):
    if v["ci_lower"] is None:
        return f"{scale * v['value']:.{digits}f}"
    return (f"{scale * v['value']:.{digits}f} {{\\scriptsize[{scale * v['ci_lower']:.{digits}f}, "
            f"{scale * v['ci_upper']:.{digits}f}]}}")


def table_acidic():
    # Organism-level Asp+Glu and Lys+Arg shares and isoelectric point per group, with organism bootstrap CIs
    r = json.load(open(RESULTS / "acidic_residue_result.json"))["groups"]
    rows = [("halophile | extreme_halophile | archaea", "Extreme halophiles", "Archaea"),
            ("halophile | moderate_halophile | archaea", "Moderate halophiles", "Archaea"),
            ("halophile | moderate_halophile | bacteria", "Moderate halophiles", "Bacteria"),
            ("halophile | slight_halophile | bacteria", "Slight halophiles", "Bacteria"),
            ("thermophile | archaea", "Thermophiles", "Archaea"),
            ("thermophile | bacteria", "Thermophiles", "Bacteria"),
            ("psychrophile | bacteria", "Psychrophiles", "Bacteria"),
            ("acidophile | bacteria", "Acidophiles", "Bacteria"),
            ("alkaliphile | bacteria", "Alkaliphiles", "Bacteria"),
            ("non-extremophile | Swiss-Prot bacteria", "Non-extremophiles", "Bacteria"),
            ("algae halophile | eukaryote", "Halophilic algae", "Eukaryotes"),
            ("algae thermophile | eukaryote", "Thermophilic algae", "Eukaryotes"),
            ("algae psychrophile | eukaryote", "Psychrophilic algae", "Eukaryotes")]
    lines = [r"\begin{tabular}{llrccc}", r"\hline",
             r"\textbf{Group} & \textbf{Domain} & \textbf{Organisms} & \textbf{Asp+Glu (\%)} & \textbf{Lys+Arg (\%)} & \textbf{pI} \\",
             r"\hline"]
    for key, name, dom in rows:
        g = r[key]
        lines.append(f"{name} & {dom} & {g['n_organisms']} & {pct_ci(g['acidic'])} & {pct_ci(g['basic'])} & "
                     f"{pct_ci(g['pI'], 2, 1)}" + r" \\")
        if key in ("halophile | slight_halophile | bacteria", "non-extremophile | Swiss-Prot bacteria"):
            lines.append(r"\hline")
    lines += [r"\hline", r"\end{tabular}"]
    path = OUT_DIR / "table_r3_acidic_residues.tex"
    path.write_text("\n".join(lines) + "\n")
    print(f"  {path}")



# ── R3: detection ───────────────────────────────────────────────────────────────────────────────────────────────────
BINARY_METRICS = [("accuracy", "Accuracy"), ("specificity", "Specificity"), ("recall", "Recall"),
                  ("f1_extremophile", "F1"), ("mcc", "MCC")]


def fig_binary():
    # A: binary models of all four representations on the clustered test set.
    # B: MCC on all test proteins vs. without the archaeal extremophiles
    bias = json.load(open(RESULTS / "binary_domain_bias_result.json"))["models"]
    fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(FULL_WIDTH, 2.9), gridspec_kw={"width_ratios": [2.2, 1]})
    x, w = np.arange(len(BINARY_METRICS)), 0.2
    for i, (key, name) in enumerate(EMBEDDERS):
        m = json.load(open(RESULTS / f"binary_{key}_clustered_result.json"))["test_metrics"]
        vals = [m[k]["value"] for k, _ in BINARY_METRICS]
        xs = x + (i - 1.5) * w
        ax_a.bar(xs, vals, w, color=style.EMBEDDER_COLORS[key], label=name)
        error_bars(ax_a, xs, vals, [m[k]["ci_lower"] for k, _ in BINARY_METRICS], [m[k]["ci_upper"] for k, _ in BINARY_METRICS])
    ax_a.set_xticks(x, [n for _, n in BINARY_METRICS], fontsize=7.5)
    ax_a.set_ylabel("Score")
    ax_a.set_ylim(0, 1)
    ax_a.legend(loc="lower center", bbox_to_anchor=(0.5, 1.02), ncol=2)
    ax_a.text(-0.14, 1.2, "A", transform=ax_a.transAxes, fontweight="bold", fontsize=10)

    xb, wb = np.arange(len(EMBEDDERS)), 0.38
    for offset, part, color, name in ((-wb / 2, "all", style.MAIN, "All test proteins"),
                                      (wb / 2, "bacteria_only", style.RAMP[0], "Bacteria only")):
        ms = [bias[key][part]["mcc"] for key, _ in EMBEDDERS]
        vals = [m["value"] for m in ms]
        ax_b.bar(xb + offset, vals, wb, color=color, label=name)
        error_bars(ax_b, xb + offset, vals, [m["ci_lower"] for m in ms], [m["ci_upper"] for m in ms])
    ax_b.set_xticks(xb, ["Compo-\nsition", "ESM2\n8M", "ESM2\n650M", "ProtT5"], fontsize=7)
    ax_b.set_ylabel("MCC")
    ax_b.set_ylim(0, 1)
    ax_b.legend(loc="lower center", bbox_to_anchor=(0.5, 1.02), ncol=1)
    ax_b.text(-0.32, 1.2, "B", transform=ax_b.transAxes, fontweight="bold", fontsize=10)
    fig.tight_layout(w_pad=1.5)
    save(fig, "fig_r7_binary")


def fig_binarization():
    # Detection on the binary test set by the binarized five-class model, the binarized six-class model and the
    # dedicated binary model (all ProtT5)
    five = json.load(open(RESULTS / "fiveclass_binarized_clustered_result.json"))["binarized_fiveclass"]
    six = json.load(open(RESULTS / "sixclass_clustered_result.json"))
    models = [(five, "Five-class, binarized", style.RAMP[0]), (six["binarized_sixclass"], "Six-class, binarized", style.RAMP[1]),
              (six["dedicated_binary"], "Binary", style.MAIN)]
    metrics = [("accuracy", "Accuracy"), ("specificity", "Specificity"), ("recall", "Recall"), ("f1", "F1"),
               ("mcc", "MCC")]
    fig, ax = plt.subplots(figsize=(FULL_WIDTH, 2.6))
    x, w = np.arange(len(metrics)), 0.26
    for i, (m, name, color) in enumerate(models):
        keys = [k if k in m else ("f1_extremophile" if k == "f1" else k) for k, _ in metrics]
        vals = [m[k]["value"] for k in keys]
        xs = x + (i - 1) * w
        ax.bar(xs, vals, w, color=color, label=name)
        error_bars(ax, xs, vals, [m[k]["ci_lower"] for k in keys], [m[k]["ci_upper"] for k in keys])
    ax.set_xticks(x, [n for _, n in metrics])
    ax.set_ylabel("Score")
    ax.set_ylim(0, 1)
    ax.legend(loc="lower center", bbox_to_anchor=(0.5, 1.02), ncol=3)
    save(fig, "fig_r8_binarization")


def fig_thermo():
    # Binary thermophile classifier (ProtT5, per protein, 30% and 20% clustering) vs. THERMO (per proteome), each next
    # to the IVYWREL baseline at the same level
    thermo = json.load(open(RESULTS / "thermo_prokaryote_result.json"))["metrics"]
    ivywrel = json.load(open(RESULTS / "ivywrel_baseline_result.json"))
    fig, ax = plt.subplots(figsize=(FULL_WIDTH, 2.6))
    labels = ["Accuracy", "Macro-F1", "F1 (thermophile)", "MCC"]
    keys_own = ["accuracy", "macro-f1_score", "f1_thermophile", "matthews-corr-coeff"]
    keys_ext = ["accuracy", "macro_f1", "f1_thermophile", "mcc"]
    x, w = np.arange(len(labels)), 0.16
    for i, (split, name, color) in enumerate([("split1_30pct", "ProtT5, 30% (per protein)", style.MAIN),
                                              ("split1_20pct", "ProtT5, 20% (per protein)", style.SECONDARY)]):
        r = json.load(open(RESULTS / f"{split}_binary_result.json"))
        vals = [r["test_metrics"][k] for k in keys_own]
        xs = x + (i - 2) * w
        ax.bar(xs, vals, w, color=color, label=name)
        error_bars(ax, xs, vals, [r["test_metrics_bootstrap"][k]["lower"] for k in keys_own],
                   [r["test_metrics_bootstrap"][k]["upper"] for k in keys_own])
    others = [(ivywrel["per_protein"]["split1_30pct"]["metrics"], "IVYWREL, 30% (per protein)", style.RAMP[0]),
              (thermo, "THERMO (per proteome)", style.BASELINE),
              (ivywrel["per_proteome"]["metrics"], "IVYWREL (per proteome)", style.light(style.BASELINE))]
    for i, (m, name, color) in enumerate(others):
        vals = [m[k]["value"] for k in keys_ext]
        xs = x + i * w
        ax.bar(xs, vals, w, color=color, label=name)
        error_bars(ax, xs, vals, [m[k]["ci_lower"] for k in keys_ext], [m[k]["ci_upper"] for k in keys_ext])
    ax.set_xticks(x, labels)
    ax.set_ylabel("Score")
    ax.set_ylim(0, 1.05)
    ax.legend(loc="lower center", bbox_to_anchor=(0.5, 1.02), ncol=3, fontsize=7)
    save(fig, "fig_r9_thermo")


# ── R4: transfer to algae ───────────────────────────────────────────────────────────────────────────────────────────
ALGAE_CLASSES = ["thermophile", "psychrophile", "halophile", "acidophile"]
ALGAE_HOMOLOGY = Path("results/homology_inference/homology_inference_result.json")


def fig_algae_recall():
    # Share of the leakage-filtered algae proteins that the binary models classify as extremophile, overall and per class
    res = json.load(open(RESULTS / "eukaryote_embedder_clustered_result.json"))["results"]
    groups = ["overall"] + ALGAE_CLASSES
    fig, ax = plt.subplots(figsize=(FULL_WIDTH, 2.6))
    x, w = np.arange(len(groups)), 0.2
    for i, (key, name) in enumerate(EMBEDDERS):
        r = res[key]
        entries = [r["overall"]] + [r["by_class"][c] for c in ALGAE_CLASSES]
        vals = [e["recall"] for e in entries]
        xs = x + (i - 1.5) * w
        ax.bar(xs, vals, w, color=style.EMBEDDER_COLORS[key], label=name)
        error_bars(ax, xs, vals, [e["ci_lower"] for e in entries], [e["ci_upper"] for e in entries])
    n = [res["prott5"]["overall"]["n_proteins"]] + [res["prott5"]["by_class"][c]["n_proteins"] for c in ALGAE_CLASSES]
    ax.set_xticks(x, [f"{g.capitalize()}\n{k:,}" for g, k in zip(["All algae"] + [label(c) for c in ALGAE_CLASSES], n)])
    ax.set_ylabel("Recall")
    ax.set_ylim(0, 1)
    ax.legend(loc="lower center", bbox_to_anchor=(0.5, 1.02), ncol=4)
    save(fig, "fig_r10_algae_recall")


def fig_algae_homology():
    # A: recall of homology-based inference (any hit to a prokaryotic extremophile protein at E <= 1e-3) vs. ProtT5.
    # B: homology recall at stricter E-value thresholds
    hom = json.load(open(ALGAE_HOMOLOGY))
    pt5 = json.load(open(RESULTS / "eukaryote_embedder_clustered_result.json"))["results"]["prott5"]
    groups = ["overall"] + ALGAE_CLASSES
    names = ["All algae"] + [label(c) for c in ALGAE_CLASSES]
    fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(FULL_WIDTH, 2.7), gridspec_kw={"width_ratios": [1.3, 1]})

    x, w = np.arange(len(groups)), 0.38
    hom_vals = [hom["results"]["0.001"][g]["recall"] for g in groups]
    pt5_entries = [pt5["overall"]] + [pt5["by_class"][c] for c in ALGAE_CLASSES]
    ax_a.bar(x - w / 2, hom_vals, w, color=style.BASELINE, label="Homology (E ≤ 0.001)")
    ax_a.bar(x + w / 2, [e["recall"] for e in pt5_entries], w, color=style.MAIN, label="ProtT5")
    error_bars(ax_a, x + w / 2, [e["recall"] for e in pt5_entries], [e["ci_lower"] for e in pt5_entries],
               [e["ci_upper"] for e in pt5_entries])
    ax_a.set_xticks(x, names, rotation=30, ha="right")
    ax_a.set_ylabel("Recall")
    ax_a.set_ylim(0, 1)
    ax_a.legend(loc="lower center", bbox_to_anchor=(0.5, 1.02), ncol=2)
    ax_a.text(-0.2, 1.1, "A", transform=ax_a.transAxes, fontweight="bold", fontsize=10)

    thresholds = [t for t in hom["thresholds"] if float(t) <= 1e-3]
    xt = np.arange(len(thresholds))
    for g, name in zip(groups, names):
        vals = [hom["results"][t][g]["recall"] for t in thresholds]
        color = style.MAIN if g == "overall" else style.CLASS_COLORS[g]
        ax_b.plot(xt, vals, marker="o", markersize=3, linewidth=1.8 if g == "overall" else 1.0, color=color, label=name)
    ax_b.set_xticks(xt, [f"$10^{{{int(np.log10(float(t)))}}}$" for t in thresholds])
    ax_b.set_xlabel("E-value threshold")
    ax_b.set_ylabel("Homology recall")
    ax_b.set_ylim(0, 0.25)
    ax_b.legend(loc="upper right", fontsize=7)
    ax_b.text(-0.25, 1.1, "B", transform=ax_b.transAxes, fontweight="bold", fontsize=10)
    fig.tight_layout(w_pad=2)
    save(fig, "fig_r11_algae_homology")


def fig_algae_confusion():
    # Five-class ProtT5 model on the filtered algae: share of each true class assigned to each predicted class
    r = json.load(open(RESULTS / "eukaryote_fiveclass_filtered_result.json"))
    classes = r["classes"]
    cm = np.array(r["eukaryote_eval"]["confusion_matrix"], dtype=float)
    rows = [classes.index(c) for c in ALGAE_CLASSES]
    cols = [classes.index(c) for c in CLASSES]
    m = cm[np.ix_(rows, cols)]
    m = m / m.sum(axis=1, keepdims=True)
    fig, ax = plt.subplots(figsize=(FULL_WIDTH * 0.62, 2.4))
    cmap = matplotlib.colors.LinearSegmentedColormap.from_list("ramp", ["#ffffff"] + style.RAMP)
    ax.imshow(m, cmap=cmap, vmin=0, vmax=0.6, aspect="auto")
    for i in range(m.shape[0]):
        for j in range(m.shape[1]):
            ax.text(j, i, f"{m[i, j]:.2f}", ha="center", va="center", fontsize=7.5,
                    color="white" if m[i, j] > 0.35 else "#1e293b",
                    fontweight="bold" if CLASSES[j] == ALGAE_CLASSES[i] else "normal")
    ax.set_xticks(range(len(CLASSES)), [label(c) for c in CLASSES], rotation=30, ha="right")
    ax.set_yticks(range(len(ALGAE_CLASSES)), [label(c) for c in ALGAE_CLASSES])
    ax.set_xlabel("Predicted class")
    ax.set_ylabel("True class")
    for spine in ax.spines.values():
        spine.set_visible(False)
    save(fig, "fig_r12_algae_confusion")


# ── Appendix ────────────────────────────────────────────────────────────────────────────────────────────────────────
def write_table(name, lines):
    path = OUT_DIR / f"{name}.tex"
    path.write_text("\n".join(lines) + "\n")
    print(f"  {path}")


def table_biotrainer():
    # Our implementation vs. biotrainer 2.1.1: accuracy, macro-F1 and MCC per task and representation
    rows = defaultdict(dict)
    with open("results/biotrainer/biotrainer_comparison.tsv") as f:
        next(f)
        for line in f:
            predictor, embedder, metric, previous, biotrainer = line.rstrip("\n").split("\t")[:5]
            rows[(predictor, embedder)][metric] = (previous, biotrainer)
    tasks = [("fiveclass_split1_30pct", "Five-class, 30\\%"), ("fiveclass_split1_20pct", "Five-class, 20\\%"),
             ("twelveclass_split2_30pct", "Twelve-class, 30\\%"), ("twelveclass_split2_20pct", "Twelve-class, 20\\%"),
             ("binary_clustered", "Binary"), ("sixclass_clustered", "Six-class")]
    names = {"prott5": "ProtT5", "esm2_650m": "ESM2 650M", "esm2_8m": "ESM2 8M", "composition": "Composition"}
    lines = [r"\begin{tabular}{llcccccc}", r"\hline",
             r" & & \multicolumn{2}{c}{\textbf{Accuracy}} & \multicolumn{2}{c}{\textbf{Macro-F1}} & \multicolumn{2}{c}{\textbf{MCC}} \\",
             r"\textbf{Task} & \textbf{Representation} & ours & biotrainer & ours & biotrainer & ours & biotrainer \\", r"\hline"]
    for key, task in tasks:
        first = True
        for emb in ("prott5", "esm2_650m", "esm2_8m", "composition"):
            if (key, emb) not in rows:
                continue
            r = rows[(key, emb)]
            cells = " & ".join(f"{r[m][0]} & {r[m][1]}" for m in ("accuracy", "macro_f1", "mcc"))
            lines.append(f"{task if first else ''} & {names[emb]} & {cells}" + r" \\")
            first = False
        lines.append(r"\hline")
    lines.append(r"\end{tabular}")
    write_table("table_a_biotrainer", lines)


def table_hierarchical():
    # Flat vs. hierarchical twelve-class model, oracle routing, and the second-step classifiers per class
    flat = multiclass_result("split2_30pct", "prott5")
    h = json.load(open(RESULTS / "split2_30pct_hierarchical_result.json"))
    metrics = [("accuracy", "Accuracy"), ("macro-f1_score", "Macro-F1"), ("matthews-corr-coeff", "MCC")]
    lines = [r"\begin{tabular}{lccc}", r"\hline",
             r"\textbf{Model} & " + " & ".join(rf"\textbf{{{n}}}" for _, n in metrics) + r" \\", r"\hline",
             "Flat twelve-class & " + " & ".join(fmt_ci(flat["test_metrics"][m], flat["test_metrics_bootstrap"][m]) for m, _ in metrics) + r" \\",
             "Hierarchical & " + " & ".join(fmt_ci(h["test_metrics"][m], h["test_metrics_bootstrap"][m]) for m, _ in metrics) + r" \\",
             "Hierarchical, oracle routing & " + " & ".join(f"{h['oracle_metrics'][m]:.3f}" for m, _ in metrics) + r" \\",
             r"\hline", r"\multicolumn{4}{l}{\textit{Second step: subtypes within a class}} \\",
             r"\textbf{Class} & \textbf{Subtypes} & \textbf{Accuracy} & \textbf{Macro-F1} \\", r"\hline"]
    for cls in ("thermophile", "halophile", "acidophile", "alkaliphile"):
        g = h["level2_sub_classifiers"][cls]
        lines.append(f"{label(cls)} & {len(g['classes'])} & {g['accuracy']:.3f} & {g['macro_f1']:.3f}" + r" \\")
    lines += [r"\hline", r"\end{tabular}"]
    write_table("table_a_hierarchical", lines)


def table_thermo_overlap():
    t = json.load(open(RESULTS / "thermo_prokaryote_result.json"))
    o, m = t["training_overlap"], t["metrics"]
    groups = [("in_thermo_training", "In THERMO's training set"), ("not_in_thermo_training", "Not in the training set"),
              ("no_bacdive_id", "No BacDive identifier")]
    lines = [r"\begin{tabular}{lrrcc}", r"\hline",
             r"\textbf{Organisms} & \textbf{Proteomes} & \textbf{Thermophiles} & \textbf{Accuracy} & \textbf{MCC} \\", r"\hline"]
    for key, name in groups:
        g = o[key]
        lines.append(f"{name} & {g['n']} & {g['n_thermophile']} & {g['accuracy']:.3f} & {g['mcc']:.3f}" + r" \\")
    lines += [r"\hline",
              f"All prokaryotes & {t['n_proteomes']} & {t['n_thermophile']} & {fmt_ci(m['accuracy']['value'], {'lower': m['accuracy']['ci_lower'], 'upper': m['accuracy']['ci_upper']})} & "
              f"{fmt_ci(m['mcc']['value'], {'lower': m['mcc']['ci_lower'], 'upper': m['mcc']['ci_upper']})}" + r" \\",
              r"\hline", r"\end{tabular}"]
    write_table("table_a_thermo_overlap", lines)


def table_literature():
    # Organisms whose growth optimum was taken directly from the primary literature (evidence type "experimental")
    rows = sorted((r for r in csv.DictReader(open(MASTER_INDEX), delimiter="\t") if r["evidence_type"] == "experimental"),
                  key=lambda r: r["taxon_name"])
    value = lambda v: f"{float(v):g}" if v else ""
    lines = [r"\begin{tabular}{lp{5cm}rrr}", r"\hline",
             r"\textbf{Organism} & \textbf{Label} & \textbf{OGT ($^\circ$C)} & \textbf{NaCl (M)} & \textbf{pH} \\",
             r"\hline"]
    for r in rows:
        name = rf"\textit{{{r['taxon_name']}}}" + ("$^{a}$" if r["proteome_status"] != "downloaded" else "")
        labels = ", ".join(r["phenotype_labels"].split("|"))
        lines.append(f"{name} & {labels} & {value(r['ogt_celsius'])} & {value(r['nacl_optimum_M'])} & "
                     f"{value(r['ph_optimum'])}" + r" \\")
    lines += [r"\hline", r"\end{tabular}"]
    write_table("table_a_literature", lines)


SUBTYPES = [("thermophile", "Thermophile"), ("hyperthermophile", "Hyperthermophile"), ("psychrophile", "Psychrophile"),
            ("slight_halophile", "Slight halophile"), ("moderate_halophile", "Moderate halophile"),
            ("extreme_halophile", "Extreme halophile"), ("acidotolerant", "Acidotolerant"),
            ("moderate_acidophile", "Moderate acidophile"), ("acidophile", "Acidophile"),
            ("alkalitolerant", "Alkalitolerant"), ("alkaliphile", "Alkaliphile"), ("extreme_alkaliphile", "Extreme alkaliphile")]


def fig_twelveclass_per_subtype():
    fig, ax = plt.subplots(figsize=(FULL_WIDTH, 2.8))
    x, w = np.arange(len(SUBTYPES)), 0.2
    for i, (key, name) in enumerate(EMBEDDERS):
        r = multiclass_result("split2_30pct", key)
        vals = [r["per_class_metrics"][s]["f1_score"] for s, _ in SUBTYPES]
        ci = [r["per_class_metrics_bootstrap"][s]["f1_score"] for s, _ in SUBTYPES]
        xs = x + (i - 1.5) * w
        ax.bar(xs, vals, w, color=style.EMBEDDER_COLORS[key], label=name)
        error_bars(ax, xs, vals, [c["lower"] for c in ci], [c["upper"] for c in ci])
    for boundary in (1.5, 2.5, 5.5, 8.5):
        ax.axvline(boundary, color="#94a3b8", linewidth=0.6, linestyle="--")
    ax.set_xticks(x, [n for _, n in SUBTYPES], rotation=40, ha="right", fontsize=7.5)
    ax.set_ylabel("F1 score")
    ax.set_ylim(0, 1)
    ax.legend(loc="lower center", bbox_to_anchor=(0.5, 1.02), ncol=4)
    save(fig, "fig_a_twelveclass_per_subtype")


def fig_crossval_per_class():
    # Per-class F1 of the five-class ProtT5 model: per-protein split vs. organism-grouped CV (mean and SD)
    fig, ax = plt.subplots(figsize=(FULL_WIDTH, 2.6))
    x, w = np.arange(len(CLASSES)), 0.26
    split = multiclass_result("split1_30pct", "prott5")
    vals = [split["per_class_metrics"][c]["f1_score"] for c in CLASSES]
    ci = [split["per_class_metrics_bootstrap"][c]["f1_score"] for c in CLASSES]
    ax.bar(x - w, vals, w, color=style.MAIN, label="Per-protein split, 30% (95% CI)")
    error_bars(ax, x - w, vals, [c["lower"] for c in ci], [c["upper"] for c in ci])
    for offset, pct, color in ((0, "30pct", style.SECONDARY), (w, "20pct", style.RAMP[0])):
        cv = json.load(open(RESULTS / f"split1_{pct}_cv5_result.json"))["cv_per_class_aggregated"]
        means = [cv[c]["f1_score"]["mean"] for c in CLASSES]
        sds = [cv[c]["f1_score"]["std"] for c in CLASSES]
        ax.bar(x + offset, means, w, color=color, label=f"Unseen organisms, {pct[:2]}% (mean ± SD)")
        ax.errorbar(x + offset, means, yerr=sds, fmt="none", ecolor="#334155", elinewidth=0.8, capsize=2)
    ax.set_xticks(x, [label(c) for c in CLASSES])
    ax.set_ylabel("F1 score")
    ax.set_ylim(0, 1)
    ax.legend(loc="lower center", bbox_to_anchor=(0.5, 1.02), ncol=3, fontsize=7)
    save(fig, "fig_a_crossval_per_class")


def fig_algae_species():
    # Recall of the binary ProtT5 model per algal species, coloured by class
    sp = json.load(open(RESULTS / "eukaryote_embedder_clustered_result.json"))["results"]["prott5"]["by_species"]
    order = {c: i for i, c in enumerate(ALGAE_CLASSES)}
    items = sorted(sp.values(), key=lambda v: (order[v["broad_class"]], -v["recall"]))
    fig, ax = plt.subplots(figsize=(FULL_WIDTH, 2.8))
    x = np.arange(len(items))
    ax.bar(x, [v["recall"] for v in items], 0.7, color=[style.CLASS_COLORS[v["broad_class"]] for v in items])
    for xi, v in zip(x, items):
        ax.text(xi, v["recall"] + 0.015, f"{v['recall']:.2f}", ha="center", fontsize=6.5)
    ax.set_xticks(x, [v["species"].replace("Chlamydomonas sp. ", "Chlamydomonas sp.\n") for v in items],
                  rotation=40, ha="right", fontsize=7, fontstyle="italic")
    ax.set_ylabel("Recall")
    ax.set_ylim(0, 1)
    handles = [plt.Rectangle((0, 0), 1, 1, color=style.CLASS_COLORS[c]) for c in ALGAE_CLASSES]
    ax.legend(handles, [label(c) for c in ALGAE_CLASSES], loc="lower center", bbox_to_anchor=(0.5, 1.02), ncol=4)
    save(fig, "fig_a_algae_species")


def fig_pfam_heatmap():
    # Top 20 Pfam domains by enrichment over the mean frequency of the five classes (same procedure as
    # scripts/interproscan_tm/analyze_interproscan.py); cells show the number of proteins with the domain
    import math
    pid2dom, desc, prot_cls = defaultdict(set), {}, {}
    with open("results/interproscan/prokaryote_results.tsv") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            pid = row["Sequence_ID"]
            prot_cls[pid] = pid.split("||")[-1]
            if row["Analysis"] == "Pfam":
                pid2dom[pid].add(row["Signature_Accession"])
                desc[row["Signature_Accession"]] = row["Signature_Description"]
    n, cnt = Counter(), defaultdict(Counter)
    for pid, cls in prot_cls.items():
        n[cls] += 1
        for d in pid2dom.get(pid, ()):
            cnt[cls][d] += 1
    eps = 0.005
    freq = lambda d, c: cnt[c][d] / n[c]
    mean = lambda d: sum(freq(d, c) for c in CLASSES) / len(CLASSES)
    enr = lambda d, c: math.log2((freq(d, c) + eps) / (mean(d) + eps))
    domains = {d for c in CLASSES for d in cnt[c]}
    top = sorted((d for d in domains if max(enr(d, c) for c in CLASSES) > 0.3),
                 key=lambda d: max(enr(d, c) for c in CLASSES), reverse=True)[:20]
    m = np.array([[enr(d, c) for c in CLASSES] for d in top])
    fig, ax = plt.subplots(figsize=(FULL_WIDTH, 5.2))
    vmax = np.abs(m).max()
    im = ax.imshow(m, cmap="RdBu_r", vmin=-vmax, vmax=vmax, aspect="auto")
    for i, d in enumerate(top):
        for j, c in enumerate(CLASSES):
            ax.text(j, i, str(cnt[c][d]), ha="center", va="center", fontsize=6.5,
                    color="white" if abs(m[i, j]) > 0.7 * vmax else "#1e293b")
    ax.set_xticks(range(len(CLASSES)), [label(c) for c in CLASSES], rotation=30, ha="right")
    ax.set_yticks(range(len(top)), [f"{d}  {desc[d][:42]}" for d in top], fontsize=6.5)
    cbar = fig.colorbar(im, ax=ax, shrink=0.5)
    cbar.set_label("log$_2$ enrichment", fontsize=8)
    save(fig, "fig_a_pfam_heatmap")

FIGURES = {
    "r1": fig_dataset_organisms,
    "r2": fig_dataset_optima,
    "data_proteome_sizes": fig_proteome_sizes,
    "r3": fig_fiveclass_per_class,
    "table_r1": table_multiclass,
    "table_a_20pct": lambda: table_multiclass("20pct", "table_a_multiclass_20pct"),
    "table_r2": table_leakage,
    "r4": fig_homology,
    "r5": fig_crossval,
    "r6": fig_archaea,
    "table_r3": table_acidic,
    "r7": fig_binary,
    "r8": fig_binarization,
    "r9": fig_thermo,
    "r10": fig_algae_recall,
    "r11": fig_algae_homology,
    "r12": fig_algae_confusion,
    "table_a_biotrainer": table_biotrainer,
    "table_a_hierarchical": table_hierarchical,
    "table_a_thermo_overlap": table_thermo_overlap,
    "table_a_literature": table_literature,
    "a_twelveclass": fig_twelveclass_per_subtype,
    "a_crossval": fig_crossval_per_class,
    "a_algae_species": fig_algae_species,
    "a_pfam": fig_pfam_heatmap,
}


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("names", nargs="*", help=f"subset of {list(FIGURES)}; default: all")
    names = parser.parse_args().names or list(FIGURES)
    setup()
    for n in names:
        print(n)
        FIGURES[n]()
