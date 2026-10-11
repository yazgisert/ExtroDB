# Binarized 5-class ProtT5 model (max probability >= 0.5 = extremophile) on the clustered binary test set, to show
# that a model without non-extremophiles cannot detect them.

import json
from collections import defaultdict

import h5py
import numpy as np
import torch
import torch.nn as nn
from sklearn.metrics import (
    accuracy_score, balanced_accuracy_score, f1_score,
    matthews_corrcoef, precision_score, recall_score, roc_auc_score,
)

SEED = 42

EXTREM_H5    = "embeddings/split1_30pct.h5"
NONEXTREM_H5 = "embeddings/nonextremophile.h5"
SPLIT_TSV    = "dataset_splits/clustered_binary/binary_split.tsv"
MODEL_PATH   = "results/training_results/fiveclass_extremophile_model.pt"
OUT_JSON     = "results/training_results/fiveclass_binarized_clustered_result.json"

BINARIZATION_THRESHOLD = 0.5
N_BOOTSTRAP = 1000


class FNN(nn.Module):
    def __init__(self, n_features=1024, n_classes=5, dropout_rate=0.25):
        super().__init__()
        self.fc1     = nn.Linear(n_features, 32)
        self.dropout = nn.Dropout(dropout_rate)
        self.fc2     = nn.Linear(32, n_classes)

    def forward(self, x):
        return self.fc2(self.dropout(torch.relu(self.fc1(x))))


def binary_metrics(y_true, y_pred):
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
    return {
        k: {"value": float(v),
            "ci_lower": float(np.percentile(boot[k], 2.5)),
            "ci_upper": float(np.percentile(boot[k], 97.5))}
        for k, v in point.items()
    }


if __name__ == "__main__":

    print("Loading the clustered binary test set ...")
    ext_ids, non_ids = [], []
    with open(SPLIT_TSV) as f:
        next(f)
        for line in f:
            protein_id, split, label = line.rstrip("\n").split("\t")
            if split == "test":
                (non_ids if label == "non_extremophile" else ext_ids).append(protein_id)
    with h5py.File(EXTREM_H5, "r") as f:
        X_ext = [f[i][:] for i in ext_ids]
    with h5py.File(NONEXTREM_H5, "r") as f:
        X_non = [f[i][:] for i in non_ids]
    print(f"  {len(X_ext):,} extremophile + {len(X_non):,} non-extremophile proteins")

    X      = np.array(X_ext + X_non, dtype=np.float32)
    y_true = np.array([1] * len(X_ext) + [0] * len(X_non))

    model = FNN()
    model.load_state_dict(torch.load(MODEL_PATH, map_location="cpu"))
    model.eval()
    with torch.no_grad():
        max_prob = torch.softmax(model(torch.from_numpy(X)), dim=1).max(dim=1).values.numpy()

    y_pred  = (max_prob >= BINARIZATION_THRESHOLD).astype(int)
    metrics = with_ci(y_true, y_pred)
    auroc   = roc_auc_score(y_true, max_prob)

    print(f"\nBinarized 5-class on prokaryote test set ({len(y_true):,} proteins, max prob >= {BINARIZATION_THRESHOLD})")
    for k, v in metrics.items():
        print(f"  {k:<18} {v['value']:.3f}  [{v['ci_lower']:.3f}, {v['ci_upper']:.3f}]")
    print(f"  {'auroc (max prob)':<18} {auroc:.3f}")

    result = {
        "task":     "fiveclass_binarized_on_prokaryote_test_set",
        "model":    "ProtT5 FNN 5-class (split1_30pct, extremophiles only), binarized by max-probability threshold",
        "binarization_threshold": BINARIZATION_THRESHOLD,
        "n_test":   int(len(y_true)),
        "test_set": "clustered binary test set (same as the binary and 6-class models)",
        "binarized_fiveclass": metrics,
        "auroc_max_prob":      float(auroc),
    }
    with open(OUT_JSON, "w") as f:
        json.dump(result, f, indent=2)
    print(f"\nResults saved to {OUT_JSON}")
