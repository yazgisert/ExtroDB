# IVYWREL baseline for thermophily (Zeldovich et al. 2007), per protein on the Split 1 test sets of train_binary.py
# (threshold fitted on the training proteins) and per proteome on the 225 proteomes of the THERMO comparison
# (published relation T_opt = 937 F - 335). Also reports the glutamate fraction of the proteomes by prediction outcome
# and class, since glutamate is one of the seven IVYWREL residues.
# Run from the repository root: python scripts/evaluation/eval_ivywrel_baseline.py

import csv
import json
from collections import defaultdict

import numpy as np
from sklearn.metrics import accuracy_score, f1_score, matthews_corrcoef

SEED = 42
N_BOOTSTRAP = 1000

IVYWREL = set("IVYWREL")
THERMOPHILE_MIN_OGT = 45.0

SPLITS = {
    "split1_30pct": ("dataset_splits/rep_split1_train.faa", "dataset_splits/rep_split1_test.faa"),
    "split1_20pct": ("dataset_splits/rep_split1_20pct_train.faa", "dataset_splits/rep_split1_20pct_test.faa"),
}
THERMO_RESULT = "results/training_results/thermo_bacdive_result.json"
MASTER_INDEX  = "searchable_resource/extremophile_master_index.tsv"
OUT_JSON      = "results/training_results/ivywrel_baseline_result.json"


def ivywrel_fraction(seq):
    return sum(aa in IVYWREL for aa in seq) / max(len(seq), 1)


def read_split_faa(path):
    # Header: >proteome||protein description||label; Split 2 labels are not used here
    fractions, labels = [], []
    with open(path) as f:
        header, seq = None, []
        for line in list(f) + [">"]:
            line = line.strip()
            if line.startswith(">"):
                if header is not None:
                    fractions.append(ivywrel_fraction("".join(seq)))
                    labels.append(1 if header.split("||")[2].strip().lower() in ("thermophile", "hyperthermophile") else 0)
                header, seq = line[1:], []
            else:
                seq.append(line)
    return np.array(fractions), np.array(labels)


def proteome_fractions(path):
    # IVYWREL and glutamate fraction over all residues of the proteome
    n_ivywrel = n_glu = n_total = 0
    with open(path) as f:
        for line in f:
            if not line.startswith(">"):
                seq = line.strip()
                n_ivywrel += sum(aa in IVYWREL for aa in seq)
                n_glu += seq.count("E")
                n_total += len(seq)
    return n_ivywrel / n_total, n_glu / n_total


def best_threshold(fractions, labels):
    candidates = np.unique(np.percentile(fractions, np.linspace(1, 99, 197)))
    scores = [matthews_corrcoef(labels, (fractions >= t).astype(int)) for t in candidates]
    return float(candidates[int(np.argmax(scores))])


def metrics(y_true, y_pred):
    return {
        "accuracy":       accuracy_score(y_true, y_pred),
        "macro_f1":       f1_score(y_true, y_pred, average="macro"),
        "f1_thermophile": f1_score(y_true, y_pred),
        "mcc":            matthews_corrcoef(y_true, y_pred),
    }


def with_ci(y_true, y_pred):
    rng = np.random.default_rng(SEED)
    boot = defaultdict(list)
    for _ in range(N_BOOTSTRAP):
        idx = rng.choice(len(y_true), size=len(y_true), replace=True)
        for k, v in metrics(y_true[idx], y_pred[idx]).items():
            boot[k].append(v)
    return {k: {"value": float(v), "ci_lower": float(np.percentile(boot[k], 2.5)),
                "ci_upper": float(np.percentile(boot[k], 97.5))} for k, v in metrics(y_true, y_pred).items()}


def report(name, result):
    print(name)
    for k, v in result.items():
        print(f"  {k:<15} {v['value']:.3f} [{v['ci_lower']:.3f}, {v['ci_upper']:.3f}]")


if __name__ == "__main__":
    result = {"per_protein": {}, "per_proteome": {}}

    for split, (faa_train, faa_test) in SPLITS.items():
        f_train, y_train = read_split_faa(faa_train)
        f_test, y_test = read_split_faa(faa_test)
        threshold = best_threshold(f_train, y_train)
        result["per_protein"][split] = {
            "threshold": threshold, "n_train": len(y_train), "n_test": len(y_test),
            "n_test_thermophile": int(y_test.sum()),
            "metrics": with_ci(y_test, (f_test >= threshold).astype(int)),
        }
        report(f"Per protein, {split} (threshold F >= {threshold:.3f} from {len(y_train):,} training proteins, "
               f"{len(y_test):,} test proteins)", result["per_protein"][split]["metrics"])

    index = {row["proteome_id"]: row
             for row in csv.DictReader(open(MASTER_INDEX), delimiter="\t") if row["domain_scope"] == "prokaryote"}
    proteomes = [p for p in json.load(open(THERMO_RESULT))["per_proteome_predictions"] if p["proteome_id"] in index]
    fractions, glutamate = np.array([proteome_fractions(index[p["proteome_id"]]["local_fasta_path"])
                                     for p in proteomes]).T
    y_true = np.array([p["true_label"] for p in proteomes])
    threshold = (THERMOPHILE_MIN_OGT + 335) / 937
    y_pred = (fractions >= threshold).astype(int)
    result["per_proteome"] = {
        "threshold": threshold, "n_proteomes": len(y_true), "n_thermophile": int(y_true.sum()),
        "confusion": {"TP": int(((y_pred == 1) & (y_true == 1)).sum()), "FN": int(((y_pred == 0) & (y_true == 1)).sum()),
                      "FP": int(((y_pred == 1) & (y_true == 0)).sum()), "TN": int(((y_pred == 0) & (y_true == 0)).sum())},
        "metrics": with_ci(y_true, y_pred),
    }
    print()
    report(f"Per proteome ({len(y_true)} proteomes, {y_true.sum()} thermophiles, F >= {threshold:.3f})",
           result["per_proteome"]["metrics"])
    print(f"  confusion       {result['per_proteome']['confusion']}")

    # Non-thermophiles are grouped by their first class; halophile|alkaliphile counts as halophile
    groups = defaultdict(list)
    for p, true, pred, glu in zip(proteomes, y_true, y_pred, glutamate):
        if true == 1:
            group = "thermophile"
        else:
            group = ("false positive " if pred == 1 else "true negative ") + \
                    index[p["proteome_id"]]["phenotype_focus"].split("|")[0]
        groups[group].append(glu)
    result["per_proteome"]["glutamate_fraction"] = {
        group: {"n": len(values), "mean": float(np.mean(values))} for group, values in sorted(groups.items())}
    print("\nGlutamate fraction of the proteomes")
    for group, values in sorted(groups.items()):
        print(f"  {group:<30} n={len(values):>3}  {np.mean(values):.3f}")

    with open(OUT_JSON, "w") as f:
        json.dump(result, f, indent=2)
    print(f"\nResults saved to {OUT_JSON}")
