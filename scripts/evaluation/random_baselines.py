# Random baseline for every test set of the thesis: each test protein gets a label drawn at random from the class
# distribution of the training set (biotrainer calls this "random sampling"). The same metrics as for the models are
# computed, with 1000 bootstrap samples for the 95% confidence intervals (seed 42). For the per-proteome comparison
# with THERMO, which has no training set of ours, the class distribution of the evaluated proteomes is used.
# All metrics are computed from the confusion matrix, which keeps the 1000 bootstrap samples fast.
# Run from the repository root (biocentral env): python scripts/evaluation/random_baselines.py

import json
from collections import Counter
from pathlib import Path

import numpy as np

SEED, N_BOOT = 42, 1000
SPLITS = Path("dataset_splits")
RESULTS = Path("results/training_results")
OUT_JSON = RESULTS / "random_baseline_result.json"


def faa_labels(path):
    return [line.rstrip("\n").split("||")[-1].strip() for line in open(path) if line.startswith(">")]


def confusion(t, p, k):
    return np.bincount(t * k + p, minlength=k * k).reshape(k, k)


def metrics_from_confusion(cm, classes, positive=None, f1_labels=None):
    n, tp = cm.sum(), np.diag(cm).astype(float)
    true, pred = cm.sum(axis=1).astype(float), cm.sum(axis=0).astype(float)
    with np.errstate(divide="ignore", invalid="ignore"):
        recall = np.nan_to_num(tp / true)
        f1 = np.nan_to_num(2 * tp / (true + pred))
        mcc = (tp.sum() * n - (pred * true).sum()) / np.sqrt((n ** 2 - (pred ** 2).sum()) * (n ** 2 - (true ** 2).sum()))
    f1_idx = [classes.index(c) for c in (f1_labels or classes)]
    m = {"accuracy": tp.sum() / n, "macro_f1": f1[f1_idx].mean(), "mcc": float(np.nan_to_num(mcc))}
    if positive is None:
        for i, c in enumerate(classes):
            m[f"f1_{c}"], m[f"recall_{c}"] = f1[i], recall[i]
    else:
        pos = classes.index(positive)
        m.update({"recall": recall[pos], "specificity": recall[1 - pos], "f1_positive": f1[pos]})
    return m


def evaluate(name, y_true, prior, positive=None, f1_labels=None, keys=None):
    classes = sorted(set(prior) | set(y_true))
    rng = np.random.default_rng(SEED)
    t = np.array([classes.index(c) for c in y_true])
    p_classes = sorted(prior)
    weights = np.array([prior[c] for c in p_classes], dtype=float)
    p = np.array([classes.index(c) for c in p_classes])[rng.choice(len(p_classes), size=len(t), p=weights / weights.sum())]
    k = len(classes)
    point = metrics_from_confusion(confusion(t, p, k), classes, positive, f1_labels)
    boot = {key: [] for key in point}
    for _ in range(N_BOOT):
        idx = rng.integers(0, len(t), len(t))
        for key, v in metrics_from_confusion(confusion(t[idx], p[idx], k), classes, positive, f1_labels).items():
            boot[key].append(v)
    keys = keys or list(point)
    result = {"n_test": len(t), "prior": {c: prior[c] / sum(prior.values()) for c in p_classes},
              "metrics": {key: {"value": float(point[key]), "ci_lower": float(np.percentile(boot[key], 2.5)),
                                "ci_upper": float(np.percentile(boot[key], 97.5))} for key in keys}}
    summary = ", ".join(f"{key} {v['value']:.3f}" for key, v in result["metrics"].items()
                        if key in ("accuracy", "macro_f1", "mcc", "recall", "specificity", "f1_positive"))
    print(f"{name:<28} n={len(t):>7,}  {summary}")
    return result


if __name__ == "__main__":
    out = {}

    # Five-class (Split 1) and twelve-class (Split 2) test sets at 30% and 20% identity
    for task, prefix in (("fiveclass_30pct", "rep_split1"), ("fiveclass_20pct", "rep_split1_20pct"),
                         ("twelveclass_30pct", "rep_split2"), ("twelveclass_20pct", "rep_split2_20pct")):
        train, test = faa_labels(SPLITS / f"{prefix}_train.faa"), faa_labels(SPLITS / f"{prefix}_test.faa")
        out[task] = evaluate(task, test, Counter(train))

    # Binary test set (extremophile vs non-extremophile), also used for the binarized five- and six-class models
    rows = [line.rstrip("\n").split("\t") for line in list(open(SPLITS / "clustered_binary/binary_split.tsv"))[1:]]
    to_binary = lambda lab: "non_extremophile" if lab == "non_extremophile" else "extremophile"
    prior = Counter(to_binary(lab) for _, s, lab in rows if s == "train")
    test = [to_binary(lab) for _, s, lab in rows if s == "test"]
    out["binary_clustered"] = evaluate("binary_clustered", test, prior, positive="extremophile")

    # Thermophile classifier per protein (thermophile vs the other four classes of Split 1)
    to_thermo = lambda lab: "thermophile" if lab == "thermophile" else "other"
    for pct, prefix in (("30pct", "rep_split1"), ("20pct", "rep_split1_20pct")):
        prior = Counter(to_thermo(lab) for lab in faa_labels(SPLITS / f"{prefix}_train.faa"))
        test = [to_thermo(lab) for lab in faa_labels(SPLITS / f"{prefix}_test.faa")]
        out[f"thermophile_{pct}"] = evaluate(f"thermophile_{pct}", test, prior, positive="thermophile")

    # Thermophile per proteome (comparison with THERMO): class distribution of the evaluated proteomes
    thermo = json.load(open(RESULTS / "thermo_prokaryote_result.json"))
    n_thermo, n_other = thermo["n_thermophile"], thermo["n_proteomes"] - thermo["n_thermophile"]
    out["thermophile_per_proteome"] = evaluate("thermophile_per_proteome", ["thermophile"] * n_thermo + ["other"] * n_other,
                                               {"thermophile": n_thermo, "other": n_other}, positive="thermophile")

    # Algae: binary models (every algal protein is an extremophile, so only the recall is defined) and the five-class
    # model (true class of the alga, predictions drawn from the five-class training distribution; macro-F1 over the
    # four classes that occur among the algae, as in eval_eukaryote_fiveclass_filtered.py)
    algae = json.load(open(RESULTS / "eukaryote_embedder_clustered_result.json"))["results"]["prott5"]
    binary_prior = Counter(to_binary(lab) for _, s, lab in rows if s == "train")
    out["algae_binary"] = {}
    for group, n in [("overall", algae["overall"]["n_proteins"])] + \
                    [(c, v["n_proteins"]) for c, v in algae["by_class"].items()]:
        out["algae_binary"][group] = evaluate(f"algae_binary {group}", ["extremophile"] * n, binary_prior,
                                              positive="extremophile", keys=["recall"])
    counts = json.load(open(RESULTS / "eukaryote_fiveclass_filtered_result.json"))["eukaryote_eval"]["class_counts"]
    present = sorted(c for c, n in counts.items() if n)
    out["algae_fiveclass"] = evaluate("algae_fiveclass", [c for c in present for _ in range(counts[c])],
                                      Counter(faa_labels(SPLITS / "rep_split1_train.faa")), f1_labels=present)

    with open(OUT_JSON, "w") as f:
        json.dump(out, f, indent=2)
    print(f"\nSaved to {OUT_JSON}")
