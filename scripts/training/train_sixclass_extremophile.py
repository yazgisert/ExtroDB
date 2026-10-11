# 6-class FNN (5 extremophile classes + non-extremophile) on the clustered 6-class split, binarized by P(non-
# extremophile) < 0.5 and compared with the dedicated binary ProtT5 model on the same test set.

import json
import os
import random
from collections import defaultdict

import h5py
import numpy as np
import torch
import torch.nn as nn
import torch.optim as optim
from sklearn.metrics import (
    accuracy_score, balanced_accuracy_score, f1_score,
    matthews_corrcoef, precision_score, recall_score,
)
from torch.utils.data import DataLoader, TensorDataset

SEED = 42
random.seed(SEED)
np.random.seed(SEED)
torch.manual_seed(SEED)

EXTREM_H5    = "embeddings/split1_30pct.h5"
NONEXTREM_H5 = "embeddings/nonextremophile.h5"
SPLIT_TSV    = "dataset_splits/clustered_binary/sixclass_split.tsv"
BINARY_MODEL = "results/training_results/binary_prott5_clustered_model.pt"
MODEL_PATH   = "results/training_results/sixclass_clustered_model.pt"
OUT_JSON     = "results/training_results/sixclass_clustered_result.json"

CLASS_NAMES = ["acidophile", "alkaliphile", "halophile", "psychrophile", "thermophile", "non_extremophile"]
LABEL2IDX   = {c: i for i, c in enumerate(CLASS_NAMES)}
NONEXTREM   = LABEL2IDX["non_extremophile"]

BINARIZATION_THRESHOLD = 0.5
N_BOOTSTRAP = 1000


class FNN(nn.Module):
    def __init__(self, n_features=1024, n_classes=6, dropout_rate=0.25):
        super().__init__()
        self.fc1     = nn.Linear(n_features, 32)
        self.dropout = nn.Dropout(dropout_rate)
        self.fc2     = nn.Linear(32, n_classes)

    def forward(self, x):
        return self.fc2(self.dropout(torch.relu(self.fc1(x))))


def load_split():
    rows = defaultdict(list)
    with open(SPLIT_TSV) as f:
        next(f)
        for line in f:
            protein_id, split, label = line.rstrip("\n").split("\t")
            rows[split].append((protein_id, LABEL2IDX[label]))

    data = {}
    with h5py.File(EXTREM_H5, "r") as ext, h5py.File(NONEXTREM_H5, "r") as non:
        for split, items in rows.items():
            X = np.array([(non if label == NONEXTREM else ext)[pid][:] for pid, label in items], dtype=np.float32)
            y = np.array([label for _, label in items])
            idx = np.random.default_rng(SEED).permutation(len(X))
            data[split] = (X[idx], y[idx])
    return data


def train_model(X_train, y_train, X_val, y_val,
                lr=1e-3, batch_size=128, max_epochs=200, patience=10, epsilon=1e-3):
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"  Device: {device}")

    X_va = torch.tensor(X_val, dtype=torch.float32).to(device)
    y_va = torch.tensor(y_val, dtype=torch.long).to(device)
    loader = DataLoader(
        TensorDataset(torch.tensor(X_train, dtype=torch.float32), torch.tensor(y_train, dtype=torch.long)),
        batch_size=batch_size, shuffle=True,
        generator=torch.Generator().manual_seed(SEED),
    )

    model     = FNN(n_classes=len(CLASS_NAMES)).to(device)
    optimizer = optim.Adam(model.parameters(), lr=lr)
    criterion = nn.CrossEntropyLoss()

    best_val_loss = float("inf")
    best_state    = None
    no_improve    = 0

    for epoch in range(max_epochs):
        model.train()
        for xb, yb in loader:
            xb, yb = xb.to(device), yb.to(device)
            optimizer.zero_grad()
            criterion(model(xb), yb).backward()
            optimizer.step()

        model.eval()
        with torch.no_grad():
            val_loss = criterion(model(X_va), y_va).item()

        if val_loss < best_val_loss - epsilon:
            best_val_loss = val_loss
            best_state    = {k: v.cpu().clone() for k, v in model.state_dict().items()}
            no_improve    = 0
        else:
            no_improve += 1
            if no_improve >= patience:
                print(f"  Early stopping at epoch {epoch + 1}")
                break

    model.load_state_dict(best_state)
    return model.to(device), device


def get_probs(model, X, device, batch_size=512):
    model.eval()
    probs = []
    with torch.no_grad():
        for i in range(0, len(X), batch_size):
            batch = torch.tensor(X[i:i + batch_size], dtype=torch.float32).to(device)
            probs.append(torch.softmax(model(batch), dim=1).cpu().numpy())
    return np.concatenate(probs)


def sixclass_metrics(y_true, y_pred):
    return {
        "accuracy":          accuracy_score(y_true, y_pred),
        "balanced_accuracy": balanced_accuracy_score(y_true, y_pred),
        "macro_f1":          f1_score(y_true, y_pred, average="macro", zero_division=0),
        "mcc":               matthews_corrcoef(y_true, y_pred),
    }


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


def with_ci(metric_fn, y_true, y_pred):
    point = metric_fn(y_true, y_pred)
    rng   = np.random.default_rng(SEED)
    boot  = defaultdict(list)
    for _ in range(N_BOOTSTRAP):
        idx = rng.choice(len(y_true), len(y_true), replace=True)
        for k, v in metric_fn(y_true[idx], y_pred[idx]).items():
            boot[k].append(v)
    return {
        k: {"value": float(v),
            "ci_lower": float(np.percentile(boot[k], 2.5)),
            "ci_upper": float(np.percentile(boot[k], 97.5))}
        for k, v in point.items()
    }


if __name__ == "__main__":

    print("Reading the clustered 6-class split ...")
    data = load_split()
    X_train, y_train = data["train"]
    X_val,   y_val   = data["val"]
    X_test,  y_test  = data["test"]
    print(f"  train: {len(y_train):,}  val: {len(y_val):,}  test: {len(y_test):,}")

    print("\nTraining 6-class model ...")
    model, device = train_model(X_train, y_train, X_val, y_val)
    torch.save(model.state_dict(), MODEL_PATH)
    print(f"  Saved to {MODEL_PATH}")

    probs   = get_probs(model, X_test, device)
    y_pred6 = probs.argmax(axis=1)

    sixclass = with_ci(sixclass_metrics, y_test, y_pred6)
    per_class_recall = {
        c: float((y_pred6[y_test == i] == i).mean()) for i, c in enumerate(CLASS_NAMES)
    }

    y_true_bin = (y_test != NONEXTREM).astype(int)
    y_pred_bin = (probs[:, NONEXTREM] < BINARIZATION_THRESHOLD).astype(int)
    binarized  = with_ci(binary_metrics, y_true_bin, y_pred_bin)

    binary_model = FNN(n_classes=2).to(device)
    binary_model.load_state_dict(torch.load(BINARY_MODEL, map_location=device))
    y_pred_binary = get_probs(binary_model, X_test, device).argmax(axis=1)
    dedicated     = with_ci(binary_metrics, y_true_bin, y_pred_binary)

    print(f"\n6-class test set ({len(y_test):,} proteins)")
    for k, v in sixclass.items():
        print(f"  {k:<18} {v['value']:.3f}")
    print("\n  Per-class recall:")
    for c, r in per_class_recall.items():
        print(f"  {c:<18} {r:.3f}")

    print(f"\nExtremophile vs. non-extremophile (threshold: P(non_extremophile) < {BINARIZATION_THRESHOLD})")
    print(f"  {'Metric':<18} {'Binarized 6-class':>18} {'Dedicated binary':>18}")
    print("  " + "-" * 56)
    for k in binarized:
        print(f"  {k:<18} {binarized[k]['value']:>18.3f} {dedicated[k]['value']:>18.3f}")

    result = {
        "task":     "sixclass_extremophile_binarized",
        "model":    "ProtT5 FNN 6-class (clustered split, all classes equal size), binarized by P(non_extremophile) threshold",
        "classes":  CLASS_NAMES,
        "binarization_threshold": BINARIZATION_THRESHOLD,
        "n_train":  int(len(X_train)),
        "n_val":    int(len(X_val)),
        "n_test":   int(len(X_test)),
        "test_set": "identical to the clustered binary test set",
        "sixclass_metrics":        sixclass,
        "sixclass_per_class_recall": per_class_recall,
        "binarized_sixclass":      binarized,
        "dedicated_binary":        dedicated,
    }
    with open(OUT_JSON, "w") as f:
        json.dump(result, f, indent=2)
    print(f"\nResults saved to {OUT_JSON}")
