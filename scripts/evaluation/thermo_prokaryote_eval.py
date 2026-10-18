# THERMO results on the prokaryotic proteomes only (the algae are evaluated separately with their own thermophile
# labels in run_thermo_eukaryote_eval.py): metrics with 1000 bootstrap CIs and the confusion matrix, from the
# per-proteome predictions of run_thermo_model.py. Also checks whether the results depend on organisms THERMO was
# trained on, by comparing BacDive strain IDs with the training strains stored in the THERMO model file.
# Run from the repository root: python scripts/evaluation/thermo_prokaryote_eval.py

import csv
import json
import pickle
import re
import warnings
from collections import defaultdict

import numpy as np
from sklearn.metrics import accuracy_score, f1_score, matthews_corrcoef, recall_score

THERMO_RESULT = "results/training_results/thermo_bacdive_result.json"
THERMO_MODEL  = "bacdive-AI/models/thermophile_data.p"
MASTER_INDEX  = "searchable_resource/extremophile_master_index.tsv"
OUT_JSON      = "results/training_results/thermo_prokaryote_result.json"

SEED = 42
N_BOOTSTRAP = 1000


def metrics(y_true, y_pred):
    return {
        "accuracy":           accuracy_score(y_true, y_pred),
        "macro_f1":           f1_score(y_true, y_pred, average="macro"),
        "f1_thermophile":     f1_score(y_true, y_pred),
        "f1_non_thermophile": f1_score(y_true, y_pred, pos_label=0),
        "mcc":                matthews_corrcoef(y_true, y_pred),
        "recall":             recall_score(y_true, y_pred),
        "specificity":        recall_score(y_true, y_pred, pos_label=0),
    }


def with_ci(y_true, y_pred):
    rng = np.random.default_rng(SEED)
    boot = defaultdict(list)
    for _ in range(N_BOOTSTRAP):
        idx = rng.integers(0, len(y_true), len(y_true))
        for k, v in metrics(y_true[idx], y_pred[idx]).items():
            boot[k].append(v)
    return {k: {"value": float(v), "ci_lower": float(np.percentile(boot[k], 2.5)),
                "ci_upper": float(np.percentile(boot[k], 97.5))} for k, v in metrics(y_true, y_pred).items()}


if __name__ == "__main__":
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")  # the model was pickled with an older scikit-learn
        training_strains = set(pickle.load(open(THERMO_MODEL, "rb"))["strains"])

    bacdive_ids, domain = {}, {}
    for row in csv.DictReader(open(MASTER_INDEX), delimiter="\t"):
        domain[row["proteome_id"]] = row["domain_scope"]
        bacdive_ids[row["proteome_id"]] = {int(i) for i in re.findall(r"bacdive\.dsmz\.de/strain/(\d+)", row["evidence_source"])}

    preds = [p for p in json.load(open(THERMO_RESULT))["per_proteome_predictions"]
             if domain.get(p["proteome_id"]) == "prokaryote"]
    y_true = np.array([p["true_label"] for p in preds])
    y_pred = np.array([p["predicted"] for p in preds])

    result = {
        "n_proteomes": len(preds),
        "n_thermophile": int(y_true.sum()),
        "confusion": {"TP": int(((y_pred == 1) & (y_true == 1)).sum()), "FN": int(((y_pred == 0) & (y_true == 1)).sum()),
                      "FP": int(((y_pred == 1) & (y_true == 0)).sum()), "TN": int(((y_pred == 0) & (y_true == 0)).sum())},
        "metrics": with_ci(y_true, y_pred),
        "training_overlap": {},
    }
    print(f"{len(preds)} prokaryotic proteomes ({y_true.sum()} thermophiles)")
    for k, v in result["metrics"].items():
        print(f"  {k:<19} {v['value']:.3f} [{v['ci_lower']:.3f}, {v['ci_upper']:.3f}]")

    groups = defaultdict(list)
    for p, t, q in zip(preds, y_true, y_pred):
        ids = bacdive_ids.get(p["proteome_id"], set())
        group = "no_bacdive_id" if not ids else "in_thermo_training" if ids & training_strains else "not_in_thermo_training"
        groups[group].append((t, q))
    print("\nTraining overlap:")
    for group in ("in_thermo_training", "not_in_thermo_training", "no_bacdive_id"):
        t, q = np.array([a for a, _ in groups[group]]), np.array([b for _, b in groups[group]])
        result["training_overlap"][group] = {"n": len(t), "n_thermophile": int(t.sum()),
                                             "accuracy": float(accuracy_score(t, q)), "mcc": float(matthews_corrcoef(t, q))}
        print(f"  {group:<24} n={len(t):>3} ({t.sum()} thermophiles)  accuracy {accuracy_score(t, q):.3f}  "
              f"MCC {matthews_corrcoef(t, q):.3f}")

    with open(OUT_JSON, "w") as f:
        json.dump(result, f, indent=2)
    print(f"\nResults saved to {OUT_JSON}")
