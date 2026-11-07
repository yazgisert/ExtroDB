# Evaluates the biotrainer models (the published checkpoints) on the same test sets as our own models and writes the
# results for the Predictor Models page, so that the page shows the values someone gets when retraining with the
# biotrainer zips. The 5- and 12-class results come straight from biotrainer's out.yml; the binary, bacteria-only,
# 6-class, algae and homology results are recomputed from the checkpoints with the evaluation functions of our scripts.
# biotrainer's FNN is the same network as ours (Linear 32, ReLU, Dropout, Linear); only the class order differs
# (biotrainer sorts the labels alphabetically), so the output rows are reordered when the weights are loaded.
# Run from the repository root (biotrainer env): python scripts/biotrainer/evaluate_biotrainer_models.py

import json
import sys
from pathlib import Path

import h5py
import numpy as np
import torch
import yaml
from safetensors.torch import load_file

SCRIPTS = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(SCRIPTS / "training"))
sys.path.insert(0, str(SCRIPTS / "evaluation"))
import eval_binary_by_domain as domain_eval
import eval_eukaryote_embedders as algae_eval
import eval_fiveclass_binarized_prokaryote as binarized_eval
import eval_homology_binary as homology_binary
import eval_homology_fiveclass as homology_fiveclass
import train_sixclass_extremophile as sixclass_eval
from train_binary_embedder import FNN, read_split, with_ci

RUNS_DIR = Path("biotrainer_runs")
OUT_JSON = Path("results/biotrainer/biotrainer_dashboard_results.json")

# biotrainer embedder name -> our embedder name, and the label shown on the website
EMBEDDERS = {"prott5": ("prott5", "ProtT5"), "esm2_650m": ("esm2_650m", "ESM2 650M"),
             "esm2_8m": ("esm2_8m", "ESM2 8M"), "composition": ("onehot", "Amino acid composition")}
MULTICLASS = ["fiveclass_split1_30pct", "fiveclass_split1_20pct", "twelveclass_split2_30pct", "twelveclass_split2_20pct"]
BINARY_CLASSES   = ["non_extremophile", "extremophile"]  # our models: 1 = extremophile
FIVECLASS_CLASSES = homology_fiveclass.CLASS_NAMES


def latest_out_yml(predictor, embedder):
    return sorted((RUNS_DIR / predictor / embedder / "output").glob("*/out.yml"), key=lambda p: p.stat().st_mtime)[-1]


def load_model(predictor, embedder, our_classes):
    out_yml = latest_out_yml(predictor, embedder)
    str2int = yaml.safe_load(open(out_yml))["derived_values"]["class_str2int"]
    weights = load_file(out_yml.parent / "hold_out_checkpoint.safetensors")
    order = [str2int[c] for c in our_classes]
    model = FNN(n_features=weights["classifier.0.weight"].shape[1], n_classes=len(our_classes))
    model.load_state_dict({"fc1.weight": weights["classifier.0.weight"], "fc1.bias": weights["classifier.0.bias"],
                           "fc2.weight": weights["classifier.3.weight"][order],
                           "fc2.bias":   weights["classifier.3.bias"][order]})
    return model.eval()


def probabilities(model, X):
    with torch.no_grad():
        return torch.softmax(model(torch.from_numpy(X)), dim=1).numpy()


def biotrainer_test_metrics(predictor, embedder):
    test = yaml.safe_load(open(latest_out_yml(predictor, embedder)))["test_results"]["test"]
    point, boot = test["inference_result"]["metrics"], {m["name"]: m for m in test["bootstrapped_metrics"]}
    # Binary runs only report f1_score (of class 0, the extremophiles)
    f1 = "macro-f1_score" if "macro-f1_score" in point else "f1_score"
    return {k: {"value": point[src], "ci_lower": boot[src]["lower"], "ci_upper": boot[src]["upper"]}
            for k, src in (("accuracy", "accuracy"), ("macro_f1", f1), ("mcc", "matthews-corr-coeff"))
            } | {"n_test": boot["accuracy"]["sample_size"]}


def pick(metrics, keys):
    return {k: metrics[k] for k in keys}


def accuracy_entry(correct):
    ci = homology_binary.accuracy_ci(correct)
    return {"value": float(correct.mean()), "ci_lower": ci[0], "ci_upper": ci[1]}


if __name__ == "__main__":
    summary = {"fiveclass": {}, "binary_clustered": {}, "binary_bacteria_only": {}, "sixclass_clustered": {},
               "algae": {}, "homology": {}}

    for predictor in MULTICLASS:
        summary["fiveclass"][predictor] = {label: biotrainer_test_metrics(predictor, bt)
                                           for bt, (_, label) in EMBEDDERS.items()}

    # Binary models on the clustered test set, complete and without the archaeal extremophiles
    domains = domain_eval.protein_info(domain_eval.proteome_domains())
    extrem, nonextrem = read_split()
    extrem_ids, nonextrem_ids = extrem["test"], nonextrem["test"]
    y_true = np.array([1] * len(extrem_ids) + [0] * len(nonextrem_ids))
    bacterial = ~np.array([domains[i][0] == "archaea" for i in extrem_ids] + [False] * len(nonextrem_ids))
    for bt, (ours, label) in EMBEDDERS.items():
        model = load_model("binary_clustered", bt, BINARY_CLASSES)
        y_pred = probabilities(model, domain_eval.features(ours, extrem_ids, nonextrem_ids)).argmax(axis=1)
        full = with_ci(y_true, y_pred)
        reported = biotrainer_test_metrics("binary_clustered", bt)
        assert abs(full["accuracy"]["value"] - reported["accuracy"]["value"]) < 1e-3, f"{bt}: checkpoint mismatch"
        summary["binary_clustered"][label] = pick(full, ("accuracy", "specificity", "recall", "f1_extremophile", "mcc"))
        summary["binary_bacteria_only"][label] = pick(with_ci(y_true[bacterial], y_pred[bacterial]),
                                                      ("accuracy", "specificity", "recall", "mcc"))
        print(f"binary {bt:<12} MCC {full['mcc']['value']:.3f} (biotrainer reports {reported['mcc']['value']:.3f})")

    # Binarized 5-class and 6-class ProtT5 models versus the dedicated binary model, on the same test set
    data = sixclass_eval.load_split()
    X_test, y_test = data["test"]
    y_true_bin = (y_test != sixclass_eval.NONEXTREM).astype(int)
    six = probabilities(load_model("sixclass_clustered", "prott5", sixclass_eval.CLASS_NAMES), X_test)
    five = probabilities(load_model("fiveclass_split1_30pct", "prott5", FIVECLASS_CLASSES), X_test)
    binary = probabilities(load_model("binary_clustered", "prott5", BINARY_CLASSES), X_test)
    keys = ("accuracy", "specificity", "recall", "mcc")
    threshold = binarized_eval.BINARIZATION_THRESHOLD
    summary["sixclass_clustered"] = {
        "Binarized 5-class": pick(binarized_eval.with_ci(y_true_bin, (five.max(axis=1) >= threshold).astype(int)), keys),
        "Binarized 6-class": pick(sixclass_eval.with_ci(sixclass_eval.binary_metrics, y_true_bin,
                                                        (six[:, sixclass_eval.NONEXTREM] < threshold).astype(int)), keys),
        "Dedicated binary":  pick(sixclass_eval.with_ci(sixclass_eval.binary_metrics, y_true_bin,
                                                        binary.argmax(axis=1)), keys),
    }

    # Algae: share of the filtered algae proteins predicted as extremophile, overall and per class
    ids, labels, _ = algae_eval.algae_proteins()
    classes = np.array([labels[i] for i in ids])
    for bt, (ours, label) in EMBEDDERS.items():
        model = load_model("binary_clustered", bt, BINARY_CLASSES)
        predicted = (probabilities(model, algae_eval.features(ours, ids))[:, 1] >= 0.5).astype(int)
        groups = {"overall": predicted, **{c: predicted[classes == c] for c in algae_eval.CLASSES}}
        recalls = {g: algae_eval.recall_with_ci(p) for g, p in groups.items()}
        summary["algae"][label] = {g: {"value": r["recall"], "ci_lower": r["ci_lower"], "ci_upper": r["ci_upper"]}
                                   for g, r in recalls.items()}
        print(f"algae  {bt:<12} recall {summary['algae'][label]['overall']['value']:.3f} ({len(ids):,} proteins)")

    # Homology-based inference (no model) against the biotrainer ProtT5 models, on all test proteins and on those
    # without a hit; the homology numbers are taken from our results, since they do not depend on the classifier
    for name, module, predictor, our_classes in (
            ("fiveclass_split1_30pct", homology_fiveclass, "fiveclass_split1_30pct", FIVECLASS_CLASSES),
            ("binary_clustered", homology_binary, "binary_clustered", BINARY_CLASSES)):
        previous = json.load(open(module.OUT_JSON))
        label, which = {}, {}
        with open(f"{module.SEARCH_DIR}/labels.tsv") as f:
            next(f)
            for line in f:
                protein_id, split, cls = line.rstrip("\n").split("\t")
                label[protein_id], which[protein_id] = cls, split
        hits = set()
        with open(f"{module.SEARCH_DIR}/hits.tsv") as f:
            for line in f:
                query, target = line.split("\t")[:2]
                if name == "binary_clustered" or (label.get(query) in our_classes and label.get(target) in our_classes):
                    hits.add(query)

        if name == "binary_clustered":
            test_ids = [i for i, s in which.items() if s == "test"]
            y = np.array([int(label[i] != "non_extremophile") for i in test_ids])
            with h5py.File(module.EXTREM_H5, "r") as ext, h5py.File(module.NONEXTREM_H5, "r") as non:
                X = np.array([(ext if y[k] else non)[i][:] for k, i in enumerate(test_ids)], dtype=np.float32)
        else:
            test_ids = [i for i, s in which.items() if s == "test" and label[i] in our_classes]
            y = np.array([our_classes.index(label[i]) for i in test_ids])
            with h5py.File(module.EXTREM_H5, "r") as f:
                X = np.array([f[i][:] for i in test_ids], dtype=np.float32)
        correct = probabilities(load_model(predictor, "prott5", our_classes), X).argmax(axis=1) == y
        has_hit = np.array([i in hits for i in test_ids])
        assert has_hit.sum() == previous["n_with_hit"], f"{name}: hit count differs from {module.OUT_JSON}"
        summary["homology"][name] = {
            "homology_accuracy":           {k: previous["homology"]["accuracy"][k] for k in ("value", "ci_lower", "ci_upper")},
            "prott5_accuracy":             accuracy_entry(correct.astype(float)),
            "prott5_accuracy_without_hit": accuracy_entry(correct[~has_hit].astype(float)),
        }
        print(f"homology {name:<22} ProtT5 {correct.mean():.3f}, without a hit {correct[~has_hit].mean():.3f}")

    OUT_JSON.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT_JSON, "w") as f:
        json.dump(summary, f, indent=2)
    print(f"\nSaved to {OUT_JSON}")
