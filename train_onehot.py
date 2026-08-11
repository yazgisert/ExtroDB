# train_onehot.py
#
# One-hot encoding (amino acid composition) baseline.
# Instead of ProtT5 embeddings, each protein is represented as a 20-dimensional
# vector of amino acid frequencies — i.e. how often each of the 20 standard
# amino acids appears in the sequence, normalized by length.
#
# Everything else is identical to train_local.py:
#   - Same FNN architecture (just n_features=20 instead of 1024)
#   - Same optimizer, early stopping, bootstrap, JSON output format
#
# This answers: how much of ProtT5's performance comes from the language model
# vs. simple sequence composition?

import json
import os
import random
from collections import defaultdict

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


# ── Reproducibility ───────────────────────────────────────────────────────────

SEED = 42
random.seed(SEED)
np.random.seed(SEED)
torch.manual_seed(SEED)

# The 20 standard amino acids (alphabetical order)
AMINO_ACIDS = list('ACDEFGHIKLMNPQRSTVWY')
AA_INDEX    = {aa: i for i, aa in enumerate(AMINO_ACIDS)}


# ── Model ────────────────────────────────────────────────────────────────────

class FNN(nn.Module):
    def __init__(self, n_features=20, n_classes=5, dropout_rate=0.25):
        super().__init__()
        self.fc1     = nn.Linear(n_features, 32)
        self.dropout = nn.Dropout(dropout_rate)
        self.fc2     = nn.Linear(32, n_classes)

    def forward(self, x):
        return self.fc2(self.dropout(torch.relu(self.fc1(x))))


# ── Data loading ──────────────────────────────────────────────────────────────

def parse_faa(filepath):
    """
    Parses a .faa file and returns {protein_id: (label, sequence)}.
    Handles multi-line sequences.
    Header format: >GCF_xxx||WP_xxx description||label
    """
    entries = {}
    current_id  = None
    current_lbl = None
    current_seq = []

    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                # Save previous entry
                if current_id is not None:
                    entries[current_id] = (current_lbl, ''.join(current_seq))
                parts      = line[1:].split('||')
                current_id  = parts[1].split()[0].strip()
                current_lbl = parts[2].strip()
                current_seq = []
            else:
                current_seq.append(line)

    # Save last entry
    if current_id is not None:
        entries[current_id] = (current_lbl, ''.join(current_seq))

    return entries


def aa_composition(sequence):
    """
    Returns a 20-dim numpy array of amino acid frequencies.
    Non-standard characters are ignored.
    """
    counts = np.zeros(20, dtype=np.float32)
    total  = 0
    for aa in sequence.upper():
        if aa in AA_INDEX:
            counts[AA_INDEX[aa]] += 1
            total += 1
    if total > 0:
        counts /= total
    return counts


def load_split(faa_train, faa_val, faa_test):
    """
    Reads the three .faa files and returns train/val/test arrays
    of amino acid composition vectors + integer labels.
    """
    train_entries = parse_faa(faa_train)
    val_entries   = parse_faa(faa_val)
    test_entries  = parse_faa(faa_test)

    all_labels = sorted(
        set(lbl for lbl, _ in train_entries.values()) |
        set(lbl for lbl, _ in val_entries.values())   |
        set(lbl for lbl, _ in test_entries.values())
    )
    label2idx = {lbl: i for i, lbl in enumerate(all_labels)}

    def build_arrays(entries):
        X, y = [], []
        for pid, (lbl, seq) in entries.items():
            X.append(aa_composition(seq))
            y.append(label2idx[lbl])
        return np.array(X, dtype=np.float32), np.array(y, dtype=np.int64)

    X_train, y_train = build_arrays(train_entries)
    X_val,   y_val   = build_arrays(val_entries)
    X_test,  y_test  = build_arrays(test_entries)

    return X_train, y_train, X_val, y_val, X_test, y_test, all_labels


# ── Training ─────────────────────────────────────────────────────────────────

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
        batch_size=batch_size, shuffle=True,
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


# ── Metrics ──────────────────────────────────────────────────────────────────

def compute_metrics(y_true, y_pred):
    corr, _ = spearmanr(y_true, y_pred)
    return {
        'accuracy':             float(accuracy_score(y_true, y_pred)),
        'balanced-accuracy':    float(balanced_accuracy_score(y_true, y_pred)),
        'macro-precision':      float(precision_score(y_true, y_pred, average='macro', zero_division=0)),
        'macro-recall':         float(recall_score(y_true, y_pred, average='macro', zero_division=0)),
        'macro-f1_score':       float(f1_score(y_true, y_pred, average='macro', zero_division=0)),
        'matthews-corr-coeff':  float(matthews_corrcoef(y_true, y_pred)),
        'spearmans-corr-coeff': float(corr) if not np.isnan(corr) else 0.0,
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
        results[name] = {
            'accuracy':  round(float(accuracy_score(y_true[mask], y_pred[mask])), 4),
            'precision': round(float(precision_score(y_true == i, y_pred == i, zero_division=0)), 4),
            'recall':    round(float(recall_score(y_true == i, y_pred == i, zero_division=0)), 4),
            'f1_score':  round(float(f1_score(y_true == i, y_pred == i, zero_division=0)), 4),
        }
    return results


def bootstrap_per_class_metrics(y_true, y_pred, class_names, n_iter=30, confidence=0.05):
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
            samples[name]['accuracy'].append(float(accuracy_score(yt[mask], yp[mask])))
            samples[name]['precision'].append(float(precision_score(yt == i, yp == i, zero_division=0)))
            samples[name]['recall'].append(float(recall_score(yt == i, yp == i, zero_division=0)))
            samples[name]['f1_score'].append(float(f1_score(yt == i, yp == i, zero_division=0)))

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


# ── Main ──────────────────────────────────────────────────────────────────────

def run_training(split_name, faa_train, faa_val, faa_test,
                 results_dir="results/training_results"):

    print(f"\n{'='*60}")
    print(f"  {split_name}  [one-hot encoding baseline]")
    print(f"{'='*60}")

    X_train, y_train, X_val, y_val, X_test, y_test, class_names = \
        load_split(faa_train, faa_val, faa_test)

    n_classes = len(class_names)
    print(f"  train: {len(X_train)}  val: {len(X_val)}  test: {len(X_test)}")
    print(f"  classes ({n_classes}): {class_names}")
    print(f"  input features: {X_train.shape[1]}  (amino acid composition)")

    print("  Training ...")
    model, history, device = train(X_train, y_train, X_val, y_val,
                                   n_classes=n_classes)

    y_pred       = predict(model, X_test, device)
    metrics      = compute_metrics(y_test, y_pred)
    boot         = bootstrap_metrics(y_test, y_pred)
    per_cls      = per_class_metrics(y_test, y_pred, class_names)
    print("  Running per-class bootstrap (n=30) ...")
    per_cls_boot = bootstrap_per_class_metrics(y_test, y_pred, class_names)

    print(f"\n  ── Test results ──")
    for k, v in metrics.items():
        print(f"    {k:28s}: {v:.4f}")
    print(f"\n  ── Per-class F1 ──")
    for cls, m in per_cls.items():
        print(f"    {cls:25s}: f1={m['f1_score']:.3f}")

    result = {
        'split':      split_name,
        'encoder':    'one_hot_composition',
        'n_features': int(X_train.shape[1]),
        'n_train':    int(len(X_train)),
        'n_val':      int(len(X_val)),
        'n_test':     int(len(X_test)),
        'classes':    class_names,
        'config': {
            'model':        'FNN',
            'architecture': 'Linear(20→32) → Dropout(0.25) → ReLU → Linear(32→n_classes)',
            'optimizer':    'Adam',
            'lr':           1e-3,
            'epsilon':      1e-3,
            'batch_size':   128,
            'max_epochs':   200,
            'patience':     10,
            'seed':         SEED,
        },
        'training_loss':               history['train_loss'],
        'validation_loss':             history['val_loss'],
        'test_metrics':                metrics,
        'test_metrics_bootstrap':      boot,
        'per_class_metrics':           per_cls,
        'per_class_metrics_bootstrap': per_cls_boot,
    }

    os.makedirs(results_dir, exist_ok=True)
    out_path = os.path.join(results_dir, f"{split_name}_onehot_result.json")
    with open(out_path, 'w') as f:
        json.dump(result, f, indent=2)
    print(f"\n  Saved → {out_path}")

    return result


if __name__ == "__main__":

    runs = [
        ("split1_30pct",
         "dataset_splits/rep_split1_train.faa",
         "dataset_splits/rep_split1_val.faa",
         "dataset_splits/rep_split1_test.faa"),

        ("split2_30pct",
         "dataset_splits/rep_split2_train.faa",
         "dataset_splits/rep_split2_val.faa",
         "dataset_splits/rep_split2_test.faa"),

        ("split1_20pct",
         "dataset_splits/rep_split1_20pct_train.faa",
         "dataset_splits/rep_split1_20pct_val.faa",
         "dataset_splits/rep_split1_20pct_test.faa"),

        ("split2_20pct",
         "dataset_splits/rep_split2_20pct_train.faa",
         "dataset_splits/rep_split2_20pct_val.faa",
         "dataset_splits/rep_split2_20pct_test.faa"),
    ]

    for split_name, faa_train, faa_val, faa_test in runs:
        run_training(split_name, faa_train, faa_val, faa_test)
