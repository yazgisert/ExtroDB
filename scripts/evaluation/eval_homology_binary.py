# Homology-based inference on the clustered binary test set (label of the best MMseqs2 hit by bit score; no hit counts
# as wrong), compared with ProtT5 overall and on proteins with and without a hit.

import json
from collections import defaultdict

import h5py
import numpy as np
import torch
import torch.nn as nn
from sklearn.metrics import (
    accuracy_score, balanced_accuracy_score, f1_score,
    matthews_corrcoef, precision_score, recall_score,
)

SEED = 42

SEARCH_DIR   = "results/test_vs_train_search_clustered"
EXTREM_H5    = "embeddings/split1_30pct.h5"
NONEXTREM_H5 = "embeddings/nonextremophile.h5"
MODEL_PATH   = "results/training_results/binary_prott5_clustered_model.pt"
OUT_JSON     = "results/training_results/homology_binary_clustered_result.json"

NO_HIT      = -1
N_BOOTSTRAP = 1000


class FNN(nn.Module):
    def __init__(self, n_features=1024, n_classes=2, dropout_rate=0.25):
        super().__init__()
        self.fc1     = nn.Linear(n_features, 32)
        self.dropout = nn.Dropout(dropout_rate)
        self.fc2     = nn.Linear(32, n_classes)

    def forward(self, x):
        return self.fc2(self.dropout(torch.relu(self.fc1(x))))


def binary_metrics(y_true, y_pred):
    # A protein without any hit gets no prediction and counts as wrong
    y_pred = np.where(y_pred == NO_HIT, 1 - y_true, y_pred)
    return {
        "accuracy":          accuracy_score(y_true, y_pred),
        "balanced_accuracy": balanced_accuracy_score(y_true, y_pred),
        "precision":         precision_score(y_true, y_pred, zero_division=0),
        "recall":            recall_score(y_true, y_pred, zero_division=0),
        "specificity":       recall_score(y_true, y_pred, pos_label=0, zero_division=0),
        "f1":                f1_score(y_true, y_pred, zero_division=0),
        "mcc":               matthews_corrcoef(y_true, y_pred),
    }


def with_ci(y_true, y_pred):
    point = binary_metrics(y_true, y_pred)
    rng   = np.random.default_rng(SEED)
    boot  = defaultdict(list)
    for _ in range(N_BOOTSTRAP):
        idx = rng.choice(len(y_true), len(y_true), replace=True)
        for k, v in binary_metrics(y_true[idx], y_pred[idx]).items():
            boot[k].append(v)
    return {k: {"value": float(v),
                "ci_lower": float(np.percentile(boot[k], 2.5)),
                "ci_upper": float(np.percentile(boot[k], 97.5))}
            for k, v in point.items()}


if __name__ == "__main__":
    label, which = {}, {}
    with open(f"{SEARCH_DIR}/labels.tsv") as f:
        next(f)
        for line in f:
            protein_id, split, cls = line.rstrip("\n").split("\t")
            label[protein_id], which[protein_id] = cls, split
    is_extremophile = {i: int(c != "non_extremophile") for i, c in label.items()}
    test_ids = [i for i, s in which.items() if s == "test"]

    best = {}
    with open(f"{SEARCH_DIR}/hits.tsv") as f:
        for line in f:
            query, target, _, _, bits = line.split("\t")[:5]
            bits = float(bits)
            if query not in best or bits > best[query][0]:
                best[query] = (bits, is_extremophile[target])

    y_true = np.array([is_extremophile[i] for i in test_ids])
    y_hbi  = np.array([best[i][1] if i in best else NO_HIT for i in test_ids])
    has_hit = y_hbi != NO_HIT

    with h5py.File(EXTREM_H5, "r") as ext, h5py.File(NONEXTREM_H5, "r") as non:
        X = np.array([(ext if is_extremophile[i] else non)[i][:] for i in test_ids], dtype=np.float32)
    model = FNN()
    model.load_state_dict(torch.load(MODEL_PATH, map_location="cpu"))
    model.eval()
    with torch.no_grad():
        y_prott5 = model(torch.from_numpy(X)).argmax(dim=1).numpy()

    result = {
        "task":       "homology-based inference vs. ProtT5 on the clustered binary test set",
        "method":     "label of the best MMseqs2 hit (bit score) among all train + val proteins; no hit counts as wrong",
        "n_test":     len(test_ids),
        "n_with_hit": int(has_hit.sum()),
        "homology":   with_ci(y_true, y_hbi),
        "prott5":     with_ci(y_true, y_prott5),
        "subsets": {
            "with_hit":    {"n": int(has_hit.sum()),
                            "homology_accuracy": float((y_hbi[has_hit] == y_true[has_hit]).mean()),
                            "prott5_accuracy":   float((y_prott5[has_hit] == y_true[has_hit]).mean())},
            "without_hit": {"n": int((~has_hit).sum()),
                            "prott5_accuracy":   float((y_prott5[~has_hit] == y_true[~has_hit]).mean())},
        },
    }

    print(f"{len(test_ids):,} test proteins, {has_hit.sum():,} ({has_hit.mean():.1%}) with a hit in train + val\n")
    print(f"{'':20}{'Homology':>10}{'ProtT5':>10}")
    for k in result["homology"]:
        print(f"  {k:<18}{result['homology'][k]['value']:>10.3f}{result['prott5'][k]['value']:>10.3f}")
    s = result["subsets"]
    print(f"\n  proteins with a hit ({s['with_hit']['n']:,}):    homology {s['with_hit']['homology_accuracy']:.3f}   ProtT5 {s['with_hit']['prott5_accuracy']:.3f}")
    print(f"  proteins without a hit ({s['without_hit']['n']:,}): homology   —     ProtT5 {s['without_hit']['prott5_accuracy']:.3f}")

    with open(OUT_JSON, "w") as f:
        json.dump(result, f, indent=2)
    print(f"\nResults saved to {OUT_JSON}")
