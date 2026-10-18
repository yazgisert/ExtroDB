# Binary thermophile (1) vs. non-thermophile (0) classifier on ProtT5 embeddings, same FNN as train_local.py.

import json
import os
import random
from collections import defaultdict

import h5py
import numpy as np
import torch
import torch.nn as nn
import torch.optim as optim
from scipy.stats import spearmanr
from sklearn.metrics import (
    accuracy_score, balanced_accuracy_score,
    precision_score, recall_score, f1_score,
    matthews_corrcoef,
)
from torch.utils.data import DataLoader, TensorDataset

SEED = 42
random.seed(SEED)
np.random.seed(SEED)
torch.manual_seed(SEED)

SPLITS = [
    ("split1_20pct",
     "embeddings/split1_20pct.h5",
     "dataset_splits/rep_split1_20pct_train.faa",
     "dataset_splits/rep_split1_20pct_val.faa",
     "dataset_splits/rep_split1_20pct_test.faa"),

    ("split2_20pct",
     "embeddings/split2_20pct.h5",
     "dataset_splits/rep_split2_20pct_train.faa",
     "dataset_splits/rep_split2_20pct_val.faa",
     "dataset_splits/rep_split2_20pct_test.faa"),

    ("split1_30pct",
     "embeddings/split1_30pct.h5",
     "dataset_splits/rep_split1_train.faa",
     "dataset_splits/rep_split1_val.faa",
     "dataset_splits/rep_split1_test.faa"),

    ("split2_30pct",
     "embeddings/split2_30pct.h5",
     "dataset_splits/rep_split2_train.faa",
     "dataset_splits/rep_split2_val.faa",
     "dataset_splits/rep_split2_test.faa"),
]

class FNN(nn.Module):
    def __init__(self, n_features=1024, n_classes=2, dropout_rate=0.25):
        super().__init__()
        self.fc1     = nn.Linear(n_features, 32)
        self.dropout = nn.Dropout(dropout_rate)
        self.fc2     = nn.Linear(32, n_classes)

    def forward(self, x):
        return self.fc2(self.dropout(torch.relu(self.fc1(x))))


def parse_faa_labels(filepath):
    labels = {}
    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                parts      = line[1:].split("||")
                protein_id = parts[1].split()[0].strip()
                label      = parts[2].strip().lower()
                # Split 2 separates thermophiles and hyperthermophiles; both count as thermophile here
                labels[protein_id] = 1 if label in ("thermophile", "hyperthermophile") else 0
    return labels


def load_split(h5_path, faa_train, faa_val, faa_test):
    train_labels = parse_faa_labels(faa_train)
    val_labels   = parse_faa_labels(faa_val)
    test_labels  = parse_faa_labels(faa_test)

    X = defaultdict(list)
    y = defaultdict(list)

    print(f"  Reading embeddings from {h5_path} ...")
    with h5py.File(h5_path, "r") as f:
        for seq_id in f:
            emb = f[seq_id][:]
            if seq_id in train_labels:
                X["train"].append(emb); y["train"].append(train_labels[seq_id])
            elif seq_id in val_labels:
                X["val"].append(emb);   y["val"].append(val_labels[seq_id])
            elif seq_id in test_labels:
                X["test"].append(emb);  y["test"].append(test_labels[seq_id])

    return (np.array(X["train"]), np.array(y["train"]),
            np.array(X["val"]),   np.array(y["val"]),
            np.array(X["test"]),  np.array(y["test"]))


def train(X_train, y_train, X_val, y_val,
          lr=1e-3, batch_size=128, max_epochs=200, patience=10, epsilon=1e-3):
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"  Device: {device}")

    X_tr = torch.tensor(X_train, dtype=torch.float32)
    y_tr = torch.tensor(y_train, dtype=torch.long)
    X_va = torch.tensor(X_val,   dtype=torch.float32).to(device)
    y_va = torch.tensor(y_val,   dtype=torch.long).to(device)

    loader = DataLoader(
        TensorDataset(X_tr, y_tr),
        batch_size=batch_size, shuffle=True,
        generator=torch.Generator().manual_seed(SEED),
    )

    model     = FNN(n_features=X_train.shape[1], n_classes=2).to(device)
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
                print(f"  Early stopping at epoch {epoch+1}")
                break

    model.load_state_dict(best_state)
    return model.to(device), device


def predict(model, X, device):
    model.eval()
    with torch.no_grad():
        logits = model(torch.tensor(X, dtype=torch.float32).to(device))
    return logits.argmax(dim=1).cpu().numpy()


def compute_metrics(y_true, y_pred):
    corr, _ = spearmanr(y_true, y_pred)
    return {
        "accuracy":             float(accuracy_score(y_true, y_pred)),
        "balanced-accuracy":    float(balanced_accuracy_score(y_true, y_pred)),
        "macro-precision":      float(precision_score(y_true, y_pred, average="macro", zero_division=0)),
        "macro-recall":         float(recall_score(y_true, y_pred, average="macro", zero_division=0)),
        "macro-f1_score":       float(f1_score(y_true, y_pred, average="macro", zero_division=0)),
        "f1_thermophile":       float(f1_score(y_true, y_pred, pos_label=1, average="binary", zero_division=0)),
        "f1_non_thermophile":   float(f1_score(y_true, y_pred, pos_label=0, average="binary", zero_division=0)),
        "matthews-corr-coeff":  float(matthews_corrcoef(y_true, y_pred)),
        "spearmans-corr-coeff": float(corr) if not np.isnan(corr) else 0.0,
    }


def bootstrap_metrics(y_true, y_pred, n_iter=30, confidence=0.05):
    rng     = np.random.default_rng(SEED)
    samples = defaultdict(list)
    n       = len(y_true)
    for _ in range(n_iter):
        idx = rng.choice(n, size=n, replace=True)
        for k, v in compute_metrics(y_true[idx], y_pred[idx]).items():
            samples[k].append(v)
    return {
        k: {
            "mean":  float(np.mean(vals)),
            "lower": float(np.percentile(vals, 100 * confidence / 2)),
            "upper": float(np.percentile(vals, 100 * (1 - confidence / 2))),
        }
        for k, vals in samples.items()
    }


def run_training(split_name, h5_path, faa_train, faa_val, faa_test,
                 results_dir="results/training_results"):

    print(f"\nBinary thermophile classifier, {split_name}")

    X_train, y_train, X_val, y_val, X_test, y_test = \
        load_split(h5_path, faa_train, faa_val, faa_test)

    n_thermo     = int(np.sum(y_train == 1))
    n_nonthermo  = int(np.sum(y_train == 0))
    print(f"  train: {len(X_train)} (thermo={n_thermo}, non-thermo={n_nonthermo})")
    print(f"  val:   {len(X_val)}  test: {len(X_test)}")

    print("  Training ...")
    model, device = train(X_train, y_train, X_val, y_val)

    y_pred = predict(model, X_test, device)
    metrics = compute_metrics(y_test, y_pred)
    boot    = bootstrap_metrics(y_test, y_pred)

    print("\n  Test results:")
    for k, v in metrics.items():
        print(f"    {k:28s}: {v:.4f}")

    result = {
        "split":      split_name,
        "task":       "binary_thermophile_vs_nonthermophile",
        "n_train":    int(len(X_train)),
        "n_val":      int(len(X_val)),
        "n_test":     int(len(X_test)),
        "classes":    ["non_thermophile", "thermophile"],
        "config": {
            "model":        "FNN",
            "architecture": "Linear(1024->32) -> Dropout(0.25) -> ReLU -> Linear(32->2)",
            "optimizer":    "Adam",
            "lr":           1e-3,
            "epsilon":      1e-3,
            "batch_size":   128,
            "max_epochs":   200,
            "patience":     10,
            "seed":         SEED,
        },
        "test_metrics":           metrics,
        "test_metrics_bootstrap": boot,
    }

    os.makedirs(results_dir, exist_ok=True)
    out_path = os.path.join(results_dir, f"{split_name}_binary_result.json")
    with open(out_path, "w") as f:
        json.dump(result, f, indent=2)
    print(f"\n  Saved to {out_path}")
    return result


if __name__ == "__main__":
    for args in SPLITS:
        run_training(*args)
