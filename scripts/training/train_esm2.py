# Trains the 5-class and 12-class classifiers on ESM-2 650M or ESM-2 8M embeddings,
# with the same FNN and training setup as train_local.py.

import argparse
import json
import os
import random

import h5py
import numpy as np
import torch
import torch.nn as nn
import torch.optim as optim
from collections import defaultdict
from scipy.stats import spearmanr
from sklearn.metrics import (
    accuracy_score, balanced_accuracy_score,
    precision_score, recall_score, f1_score,
    matthews_corrcoef,
)
from torch.utils.data import DataLoader, TensorDataset


SEED = 42
N_BOOTSTRAP = 1000  # was 30 before 2026-10-18; the models are unchanged
random.seed(SEED)
np.random.seed(SEED)
torch.manual_seed(SEED)


class FNN(nn.Module):
    def __init__(self, n_features=1280, n_classes=5, dropout_rate=0.25):
        super().__init__()
        self.fc1     = nn.Linear(n_features, 32)
        self.dropout = nn.Dropout(dropout_rate)
        self.fc2     = nn.Linear(32, n_classes)

    def forward(self, x):
        return self.fc2(self.dropout(torch.relu(self.fc1(x))))


def parse_faa_labels(filepath):
    # Header format: >GCF_xxx||WP_xxx description||label
    labels = {}
    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                parts = line[1:].split('||')
                protein_id = parts[1].split()[0].strip()
                label      = parts[2].strip()
                labels[protein_id] = label
    return labels


def load_split(h5_path, faa_train, faa_val, faa_test):
    train_labels = parse_faa_labels(faa_train)
    val_labels   = parse_faa_labels(faa_val)
    test_labels  = parse_faa_labels(faa_test)

    all_labels  = sorted(set(train_labels.values()) |
                         set(val_labels.values())   |
                         set(test_labels.values()))
    label2idx   = {lbl: i for i, lbl in enumerate(all_labels)}

    X = defaultdict(list)
    y = defaultdict(list)

    print(f"  Reading embeddings from {h5_path} ...")
    with h5py.File(h5_path, 'r') as f:
        for seq_id in f:
            emb = f[seq_id][:]
            if seq_id in train_labels:
                X['train'].append(emb);  y['train'].append(label2idx[train_labels[seq_id]])
            elif seq_id in val_labels:
                X['val'].append(emb);    y['val'].append(label2idx[val_labels[seq_id]])
            elif seq_id in test_labels:
                X['test'].append(emb);   y['test'].append(label2idx[test_labels[seq_id]])

    return (np.array(X['train']), np.array(y['train']),
            np.array(X['val']),   np.array(y['val']),
            np.array(X['test']),  np.array(y['test']),
            all_labels)


def train(X_train, y_train, X_val, y_val, n_classes,
          lr=1e-3, batch_size=128, max_epochs=200, patience=10, epsilon=1e-3):
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    print(f"  Device: {device}")

    X_tr = torch.tensor(X_train, dtype=torch.float32)
    y_tr = torch.tensor(y_train, dtype=torch.long)
    X_va = torch.tensor(X_val,   dtype=torch.float32).to(device)
    y_va = torch.tensor(y_val,   dtype=torch.long).to(device)

    loader = DataLoader(
        TensorDataset(X_tr, y_tr),
        batch_size=batch_size,
        shuffle=True,
        generator=torch.Generator().manual_seed(SEED),
    )

    model     = FNN(n_features=X_train.shape[1], n_classes=n_classes).to(device)
    optimizer = optim.Adam(model.parameters(), lr=lr)
    criterion = nn.CrossEntropyLoss()

    history       = {'train_loss': [], 'val_loss': []}
    best_val_loss = float('inf')
    best_state    = None
    no_improve    = 0

    for epoch in range(max_epochs):
        model.train()
        epoch_loss = 0.0
        for xb, yb in loader:
            xb, yb = xb.to(device), yb.to(device)
            optimizer.zero_grad()
            loss = criterion(model(xb), yb)
            loss.backward()
            optimizer.step()
            epoch_loss += loss.item() * len(xb)
        train_loss = epoch_loss / len(X_train)

        model.eval()
        with torch.no_grad():
            val_loss = criterion(model(X_va), y_va).item()

        history['train_loss'].append(round(train_loss, 6))
        history['val_loss'].append(round(val_loss, 6))

        if (epoch + 1) % 10 == 0:
            print(f"  Epoch {epoch+1:3d}  train_loss={train_loss:.4f}  val_loss={val_loss:.4f}")

        if val_loss < best_val_loss - epsilon:
            best_val_loss = val_loss
            best_state    = {k: v.cpu().clone() for k, v in model.state_dict().items()}
            no_improve    = 0
        else:
            no_improve += 1
            if no_improve >= patience:
                print(f"  Early stopping at epoch {epoch+1}  (best val loss: {best_val_loss:.4f})")
                break

    model.load_state_dict(best_state)
    return model.to(device), history, device


def predict(model, X, device):
    model.eval()
    with torch.no_grad():
        logits = model(torch.tensor(X, dtype=torch.float32).to(device))
    return logits.argmax(dim=1).cpu().numpy()


def compute_metrics(y_true, y_pred):
    corr, _ = spearmanr(y_true, y_pred)
    return {
        'accuracy':              float(accuracy_score(y_true, y_pred)),
        'balanced-accuracy':     float(balanced_accuracy_score(y_true, y_pred)),
        'macro-precision':       float(precision_score(y_true, y_pred, average='macro', zero_division=0)),
        'macro-recall':          float(recall_score(y_true, y_pred, average='macro', zero_division=0)),
        'macro-f1_score':        float(f1_score(y_true, y_pred, average='macro', zero_division=0)),
        'matthews-corr-coeff':   float(matthews_corrcoef(y_true, y_pred)),
        'spearmans-corr-coeff':  float(corr) if not np.isnan(corr) else 0.0,
    }


def bootstrap_metrics(y_true, y_pred, n_iter=N_BOOTSTRAP, confidence=0.05):
    rng     = np.random.default_rng(SEED)
    samples = defaultdict(list)
    n       = len(y_true)

    for _ in range(n_iter):
        idx = rng.choice(n, size=n, replace=True)
        for k, v in compute_metrics(y_true[idx], y_pred[idx]).items():
            samples[k].append(v)

    return {
        k: {
            'mean':  float(np.mean(vals)),
            'lower': float(np.percentile(vals, 100 * confidence / 2)),
            'upper': float(np.percentile(vals, 100 * (1 - confidence / 2))),
        }
        for k, vals in samples.items()
    }


def per_class_metrics(y_true, y_pred, class_names):
    results = {}
    for i, name in enumerate(class_names):
        mask = y_true == i
        if mask.sum() == 0:
            continue
        acc  = accuracy_score(y_true[mask], y_pred[mask])
        prec = precision_score(y_true == i, y_pred == i, zero_division=0)
        rec  = recall_score(y_true == i, y_pred == i, zero_division=0)
        f1   = f1_score(y_true == i, y_pred == i, zero_division=0)
        results[name] = {
            'accuracy':  round(float(acc),  4),
            'precision': round(float(prec), 4),
            'recall':    round(float(rec),  4),
            'f1_score':  round(float(f1),   4),
        }
    return results


def bootstrap_per_class_metrics(y_true, y_pred, class_names, n_iter=N_BOOTSTRAP, confidence=0.05):
    rng     = np.random.default_rng(SEED)
    n       = len(y_true)
    samples = {name: {'accuracy': [], 'precision': [], 'recall': [], 'f1_score': []}
               for name in class_names}

    for _ in range(n_iter):
        idx = rng.choice(n, size=n, replace=True)
        yt, yp = y_true[idx], y_pred[idx]
        for i, name in enumerate(class_names):
            mask = yt == i
            if mask.sum() == 0:
                continue
            samples[name]['accuracy'].append(
                float(accuracy_score(yt[mask], yp[mask])))
            samples[name]['precision'].append(
                float(precision_score(yt == i, yp == i, zero_division=0)))
            samples[name]['recall'].append(
                float(recall_score(yt == i, yp == i, zero_division=0)))
            samples[name]['f1_score'].append(
                float(f1_score(yt == i, yp == i, zero_division=0)))

    def ci(vals):
        return {
            'mean':  float(np.mean(vals)),
            'lower': float(np.percentile(vals, 100 * confidence / 2)),
            'upper': float(np.percentile(vals, 100 * (1 - confidence / 2))),
        }

    return {
        name: {metric: ci(vals) for metric, vals in metrics.items() if vals}
        for name, metrics in samples.items()
    }


def run_training(split_name, h5_path, faa_train, faa_val, faa_test, suffix,
                 results_dir="results/training_results"):

    print(f"\n{split_name}")

    X_train, y_train, X_val, y_val, X_test, y_test, class_names = \
        load_split(h5_path, faa_train, faa_val, faa_test)

    n_classes = len(class_names)
    print(f"  train: {len(X_train)}  val: {len(X_val)}  test: {len(X_test)}")
    print(f"  classes ({n_classes}): {class_names}")

    print("Training ...")
    model, history, device = train(X_train, y_train, X_val, y_val,
                                   n_classes=n_classes)

    y_pred        = predict(model, X_test, device)
    metrics       = compute_metrics(y_test, y_pred)
    boot          = bootstrap_metrics(y_test, y_pred)
    per_cls       = per_class_metrics(y_test, y_pred, class_names)
    per_cls_boot  = bootstrap_per_class_metrics(y_test, y_pred, class_names)

    print("\n  Test results:")
    for k, v in metrics.items():
        print(f"    {k:28s}: {v:.4f}")
    print("\n  Per-class accuracy:")
    for cls, m in per_cls.items():
        print(f"    {cls:20s}: acc={m['accuracy']:.3f}  f1={m['f1_score']:.3f}")

    result = {
        'split':      split_name,
        'n_train':    int(len(X_train)),
        'n_val':      int(len(X_val)),
        'n_test':     int(len(X_test)),
        'classes':    class_names,
        'config': {
            'model':        'FNN',
            'architecture': f'Linear({X_train.shape[1]}→32) → Dropout(0.25) → ReLU → Linear(32→n_classes)',
            'optimizer':    'Adam',
            'lr':           1e-3,
            'epsilon':      1e-3,
            'batch_size':   128,
            'max_epochs':   200,
            'patience':     10,
            'seed':         SEED,
        },
        'training_loss':          history['train_loss'],
        'validation_loss':        history['val_loss'],
        'test_metrics':                metrics,
        'test_metrics_bootstrap':      boot,
        'per_class_metrics':           per_cls,
        'per_class_metrics_bootstrap': per_cls_boot,
    }

    os.makedirs(results_dir, exist_ok=True)
    out_path = os.path.join(results_dir, f"{split_name}{suffix}_result.json")
    with open(out_path, 'w') as f:
        json.dump(result, f, indent=2)
    print(f"\n  Saved to {out_path}")

    return result


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--embedder", choices=["esm2_650m", "esm2_8m"], default="esm2_650m")
    parser.add_argument("--results-dir", default="results/training_results")
    args = parser.parse_args()
    suffix = {"esm2_650m": "_esm2", "esm2_8m": "_esm2_8m"}[args.embedder]

    runs = [
        ("split1_30pct",
         f"embeddings/split1_30pct{suffix}.h5",
         "dataset_splits/rep_split1_train.faa",
         "dataset_splits/rep_split1_val.faa",
         "dataset_splits/rep_split1_test.faa"),

        ("split2_30pct",
         f"embeddings/split2_30pct{suffix}.h5",
         "dataset_splits/rep_split2_train.faa",
         "dataset_splits/rep_split2_val.faa",
         "dataset_splits/rep_split2_test.faa"),

        ("split1_20pct",
         f"embeddings/split1_20pct{suffix}.h5",
         "dataset_splits/rep_split1_20pct_train.faa",
         "dataset_splits/rep_split1_20pct_val.faa",
         "dataset_splits/rep_split1_20pct_test.faa"),

        ("split2_20pct",
         f"embeddings/split2_20pct{suffix}.h5",
         "dataset_splits/rep_split2_20pct_train.faa",
         "dataset_splits/rep_split2_20pct_val.faa",
         "dataset_splits/rep_split2_20pct_test.faa"),
    ]

    for split_name, h5_path, faa_train, faa_val, faa_test in runs:
        run_training(split_name, h5_path, faa_train, faa_val, faa_test, suffix, args.results_dir)
