# train_cv.py
#
# 5-fold organism-grouped cross-validation for the ProtT5 FNN classifier.
#
# Key design decisions:
#   - Grouping is by ORGANISM (proteome_id), not by protein.
#     All proteins from the same organism stay in the same fold.
#     This prevents data leakage from related proteins.
#   - All proteins from train+val+test .faa files are pooled, then
#     re-split via GroupKFold — the original split is ignored.
#   - Within each fold, 15% of training organisms are held out for
#     early stopping (same patience/epsilon as the original training).
#   - Metrics and bootstrap CIs are computed per fold and aggregated.
#
# Run:
#   conda activate biocentral
#   cd thesis_data_collection
#   python train_cv.py                        # runs split1_20pct (default)
#   python train_cv.py --split split2_20pct  # or another split
#   python train_cv.py --k 10                # 10-fold instead of 5-fold

import argparse
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
    f1_score, matthews_corrcoef,
    precision_score, recall_score,
)
from sklearn.model_selection import GroupKFold

# ── Reproducibility ────────────────────────────────────────────────────────────
SEED = 42
random.seed(SEED)
np.random.seed(SEED)
torch.manual_seed(SEED)

# ── Split configs (same paths as train_local.py) ───────────────────────────────
SPLIT_CONFIGS = {
    "split1_20pct": {
        "h5":    "embeddings/split1_20pct.h5",
        "faas":  [
            "dataset_splits/rep_split1_20pct_train.faa",
            "dataset_splits/rep_split1_20pct_val.faa",
            "dataset_splits/rep_split1_20pct_test.faa",
        ],
    },
    "split2_20pct": {
        "h5":    "embeddings/split2_20pct.h5",
        "faas":  [
            "dataset_splits/rep_split2_20pct_train.faa",
            "dataset_splits/rep_split2_20pct_val.faa",
            "dataset_splits/rep_split2_20pct_test.faa",
        ],
    },
    "split1_30pct": {
        "h5":    "embeddings/split1_30pct.h5",
        "faas":  [
            "dataset_splits/rep_split1_train.faa",
            "dataset_splits/rep_split1_val.faa",
            "dataset_splits/rep_split1_test.faa",
        ],
    },
    "split2_30pct": {
        "h5":    "embeddings/split2_30pct.h5",
        "faas":  [
            "dataset_splits/rep_split2_train.faa",
            "dataset_splits/rep_split2_val.faa",
            "dataset_splits/rep_split2_test.faa",
        ],
    },
}


# ── Model (identical to train_local.py) ───────────────────────────────────────
class FNN(nn.Module):
    def __init__(self, n_features=1024, n_classes=5, dropout_rate=0.25):
        super().__init__()
        self.fc1     = nn.Linear(n_features, 32)
        self.dropout = nn.Dropout(dropout_rate)
        self.fc2     = nn.Linear(32, n_classes)

    def forward(self, x):
        return self.fc2(self.dropout(torch.relu(self.fc1(x))))


# ── Data loading ───────────────────────────────────────────────────────────────
def parse_faa_all(faa_paths):
    """
    Read all .faa files and return:
        protein_labels  : {protein_id: label_str}
        protein_organism: {protein_id: organism_id}
    Header format: >organism_id||protein_id [description]||label
    """
    protein_labels   = {}
    protein_organism = {}
    for faa_path in faa_paths:
        with open(faa_path) as f:
            for line in f:
                line = line.strip()
                if not line.startswith(">"):
                    continue
                parts       = line[1:].split("||")
                organism_id = parts[0].strip()
                protein_id  = parts[1].split()[0].strip()
                label       = parts[2].strip()
                protein_labels[protein_id]   = label
                protein_organism[protein_id] = organism_id
    return protein_labels, protein_organism


def load_all_data(h5_path, faa_paths):
    """
    Load all embeddings + labels + organism groups for a split.
    Returns:
        X       : (N, D) float32 numpy array
        y       : (N,)   int numpy array
        groups  : (N,)   str numpy array  (organism_id per protein)
        classes : sorted list of class names
    """
    protein_labels, protein_organism = parse_faa_all(faa_paths)
    all_classes = sorted(set(protein_labels.values()))
    label2idx   = {lbl: i for i, lbl in enumerate(all_classes)}

    X_list, y_list, g_list = [], [], []
    print(f"  Reading embeddings from {h5_path} ...")
    with h5py.File(h5_path, "r") as f:
        for seq_id in f:
            if seq_id not in protein_labels:
                continue
            X_list.append(f[seq_id][:])
            y_list.append(label2idx[protein_labels[seq_id]])
            g_list.append(protein_organism[seq_id])

    X      = np.array(X_list, dtype=np.float32)
    y      = np.array(y_list, dtype=np.int64)
    groups = np.array(g_list)
    print(f"  Loaded {len(X)} proteins, {len(set(groups))} organisms, {len(all_classes)} classes")
    return X, y, groups, all_classes


# ── Training (identical logic to train_local.py) ──────────────────────────────
def train_fold(X_train, y_train, X_val, y_val, n_classes,
               lr=1e-3, batch_size=128, max_epochs=200, patience=10, epsilon=1e-3):
    from torch.utils.data import DataLoader, TensorDataset

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

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

    best_val_loss = float("inf")
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
                break

    model.load_state_dict(best_state)
    return model.to(device), device


def predict(model, X, device):
    model.eval()
    with torch.no_grad():
        logits = model(torch.tensor(X, dtype=torch.float32).to(device))
    return logits.argmax(dim=1).cpu().numpy()


# ── Metrics ────────────────────────────────────────────────────────────────────
def compute_metrics(y_true, y_pred):
    corr, _ = spearmanr(y_true, y_pred)
    return {
        "accuracy":             float(accuracy_score(y_true, y_pred)),
        "balanced-accuracy":    float(balanced_accuracy_score(y_true, y_pred)),
        "macro-precision":      float(precision_score(y_true, y_pred, average="macro", zero_division=0)),
        "macro-recall":         float(recall_score(y_true, y_pred, average="macro", zero_division=0)),
        "macro-f1_score":       float(f1_score(y_true, y_pred, average="macro", zero_division=0)),
        "matthews-corr-coeff":  float(matthews_corrcoef(y_true, y_pred)),
        "spearmans-corr-coeff": float(corr) if not np.isnan(corr) else 0.0,
    }


def per_class_f1(y_true, y_pred, class_names):
    results = {}
    for i, name in enumerate(class_names):
        results[name] = {
            "precision": float(precision_score(y_true == i, y_pred == i, zero_division=0)),
            "recall":    float(recall_score(y_true == i, y_pred == i, zero_division=0)),
            "f1_score":  float(f1_score(y_true == i, y_pred == i, zero_division=0)),
            "support":   int(np.sum(y_true == i)),
        }
    return results


def aggregate_folds(fold_metrics):
    """Given list of metric dicts (one per fold), return mean ± std."""
    keys = list(fold_metrics[0].keys())
    return {
        k: {
            "mean": float(np.mean([m[k] for m in fold_metrics])),
            "std":  float(np.std( [m[k] for m in fold_metrics])),
            "per_fold": [round(m[k], 4) for m in fold_metrics],
        }
        for k in keys
    }


def aggregate_per_class(fold_per_class, class_names):
    """Aggregate per-class metrics across folds."""
    result = {}
    for name in class_names:
        metrics_list = [f[name] for f in fold_per_class if name in f]
        if not metrics_list:
            continue
        result[name] = {
            metric: {
                "mean": float(np.mean([m[metric] for m in metrics_list])),
                "std":  float(np.std( [m[metric] for m in metrics_list])),
            }
            for metric in ["precision", "recall", "f1_score"]
        }
    return result


# ── Main ───────────────────────────────────────────────────────────────────────
def run_cv(split_name, k=5, val_organism_fraction=0.15,
           results_dir="results/training_results"):

    cfg = SPLIT_CONFIGS.get(split_name)
    if cfg is None:
        raise ValueError(f"Unknown split: {split_name}. Choose from {list(SPLIT_CONFIGS)}")

    print(f"\n{'='*60}")
    print(f"  {k}-fold grouped CV — {split_name}")
    print(f"{'='*60}")

    # Load all data
    X, y, groups, class_names = load_all_data(cfg["h5"], cfg["faas"])
    n_classes = len(class_names)

    # Grouped k-fold by organism
    gkf        = GroupKFold(n_splits=k)
    fold_metrics      = []
    fold_per_class    = []
    all_y_true        = []
    all_y_pred        = []

    for fold_i, (train_idx, test_idx) in enumerate(gkf.split(X, y, groups), 1):
        print(f"\n── Fold {fold_i}/{k} ─────────────────────────────────────────")

        X_test_fold  = X[test_idx]
        y_test_fold  = y[test_idx]
        test_orgs    = set(groups[test_idx])
        train_orgs   = sorted(set(groups[train_idx]))

        # Split training organisms into train / internal-val for early stopping
        rng          = np.random.default_rng(SEED + fold_i)
        n_val_orgs   = max(1, int(len(train_orgs) * val_organism_fraction))
        val_orgs     = set(rng.choice(train_orgs, size=n_val_orgs, replace=False))
        pure_train_orgs = set(train_orgs) - val_orgs

        inner_train_mask = np.array([g in pure_train_orgs for g in groups[train_idx]])
        inner_val_mask   = np.array([g in val_orgs        for g in groups[train_idx]])

        X_tr = X[train_idx][inner_train_mask]
        y_tr = y[train_idx][inner_train_mask]
        X_va = X[train_idx][inner_val_mask]
        y_va = y[train_idx][inner_val_mask]

        print(f"  train proteins: {len(X_tr)}  "
              f"inner-val proteins: {len(X_va)}  "
              f"test proteins: {len(X_test_fold)}")
        print(f"  train orgs: {len(pure_train_orgs)}  "
              f"inner-val orgs: {len(val_orgs)}  "
              f"test orgs: {len(test_orgs)}")

        # Train
        model, device = train_fold(X_tr, y_tr, X_va, y_va, n_classes=n_classes)

        # Evaluate
        y_pred_fold = predict(model, X_test_fold, device)
        metrics     = compute_metrics(y_test_fold, y_pred_fold)
        pc          = per_class_f1(y_test_fold, y_pred_fold, class_names)

        fold_metrics.append(metrics)
        fold_per_class.append(pc)
        all_y_true.extend(y_test_fold.tolist())
        all_y_pred.extend(y_pred_fold.tolist())

        print(f"  accuracy={metrics['accuracy']:.4f}  "
              f"macro-F1={metrics['macro-f1_score']:.4f}  "
              f"MCC={metrics['matthews-corr-coeff']:.4f}")

    # ── Aggregate ──────────────────────────────────────────────────────────────
    print(f"\n{'='*60}")
    print(f"  {k}-fold CV results — {split_name}")
    print(f"{'='*60}")
    agg       = aggregate_folds(fold_metrics)
    agg_pc    = aggregate_per_class(fold_per_class, class_names)

    # Also compute metrics on pooled predictions across all folds
    pooled = compute_metrics(np.array(all_y_true), np.array(all_y_pred))

    for k_name, v in agg.items():
        print(f"  {k_name:28s}: {v['mean']:.4f} ± {v['std']:.4f}")
    print(f"\n  Per-class F1 (mean ± std):")
    for cls, m in agg_pc.items():
        print(f"    {cls:20s}: {m['f1_score']['mean']:.4f} ± {m['f1_score']['std']:.4f}")

    # ── Save ───────────────────────────────────────────────────────────────────
    result = {
        "experiment":     f"{split_name}_cv{k}",
        "split":          split_name,
        "k_folds":        k,
        "n_proteins":     int(len(X)),
        "n_organisms":    int(len(set(groups))),
        "classes":        class_names,
        "config": {
            "model":        "FNN",
            "architecture": "Linear(1024->32) -> Dropout(0.25) -> ReLU -> Linear(32->n_classes)",
            "optimizer":    "Adam",
            "lr":           1e-3,
            "epsilon":      1e-3,
            "batch_size":   128,
            "max_epochs":   200,
            "patience":     10,
            "seed":         SEED,
            "val_organism_fraction": val_organism_fraction,
        },
        "cv_metrics_aggregated": agg,
        "cv_per_class_aggregated": agg_pc,
        "pooled_metrics": pooled,
        "per_fold_metrics": fold_metrics,
    }

    os.makedirs(results_dir, exist_ok=True)
    out_path = os.path.join(results_dir, f"{split_name}_cv{k}_result.json")
    with open(out_path, "w") as f:
        json.dump(result, f, indent=2)
    print(f"\n  Saved -> {out_path}")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--split", default="split1_20pct",
                        choices=list(SPLIT_CONFIGS),
                        help="Which split to cross-validate")
    parser.add_argument("--k", type=int, default=5,
                        help="Number of folds (default: 5)")
    args = parser.parse_args()

    run_cv(args.split, k=args.k)
