# train_hierarchical.py
#
# Hierarchical predictor for fine-grained extremophile classification.
#
# Architecture:
#   Level 1 (coarse): 5-class FNN trained on split1_30pct ProtT5 embeddings
#                     → acidophile | alkaliphile | halophile | psychrophile | thermophile
#
#   Level 2 (fine):   4 sub-classifiers trained on split2_30pct ProtT5 embeddings,
#                     one per group with multiple sub-classes:
#                     → halophile   : slight / moderate / extreme halophile
#                     → acidophile  : acidotolerant / moderate_acidophile / acidophile
#                     → alkaliphile : alkalitolerant / alkaliphile / extreme_alkaliphile
#                     → thermophile : thermophile / hyperthermophile
#                     (psychrophile has no sub-classes — Level 1 is the final answer)
#
# At inference: Level 1 routes each sample to its group → Level 2 predicts within the group.
#
# Evaluation: hierarchical predictions on split2_30pct test set (12-class),
#             compared to the flat 12-class baseline. Also reports oracle accuracy
#             (what we'd achieve if Level 1 routing were perfect).
#
# Run:
#   conda activate biocentral
#   cd thesis_data_collection
#   python train_hierarchical.py

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


# ── Reproducibility ───────────────────────────────────────────────────────────

SEED = 42
random.seed(SEED)
np.random.seed(SEED)
torch.manual_seed(SEED)


# ── Hierarchy definition ──────────────────────────────────────────────────────

# Fine-grained label → coarse group
COARSE_MAP = {
    'acidophile':          'acidophile',
    'acidotolerant':       'acidophile',
    'moderate_acidophile': 'acidophile',
    'alkaliphile':         'alkaliphile',
    'alkalitolerant':      'alkaliphile',
    'extreme_alkaliphile': 'alkaliphile',
    'extreme_halophile':   'halophile',
    'moderate_halophile':  'halophile',
    'slight_halophile':    'halophile',
    'thermophile':         'thermophile',
    'hyperthermophile':    'thermophile',
    'psychrophile':        'psychrophile',
}

# All sub-classes per coarse group (sorted for reproducibility)
SUBGROUPS = {
    'acidophile':  sorted(['acidotolerant', 'moderate_acidophile', 'acidophile']),
    'alkaliphile': sorted(['alkalitolerant', 'alkaliphile', 'extreme_alkaliphile']),
    'halophile':   sorted(['slight_halophile', 'moderate_halophile', 'extreme_halophile']),
    'thermophile': sorted(['thermophile', 'hyperthermophile']),
    'psychrophile': ['psychrophile'],   # leaf — no sub-classifier needed
}

COARSE_CLASSES = sorted(SUBGROUPS.keys())
FINE_CLASSES   = sorted(COARSE_MAP.keys())   # all 12 fine-grained labels


# ── Model ─────────────────────────────────────────────────────────────────────

class FNN(nn.Module):
    def __init__(self, n_features=1024, n_classes=5, dropout_rate=0.25):
        super().__init__()
        self.fc1     = nn.Linear(n_features, 32)
        self.dropout = nn.Dropout(dropout_rate)
        self.fc2     = nn.Linear(32, n_classes)

    def forward(self, x):
        return self.fc2(self.dropout(torch.relu(self.fc1(x))))


# ── Data loading ──────────────────────────────────────────────────────────────

def parse_faa_labels(filepath):
    """Returns {protein_id: label}."""
    labels = {}
    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                parts = line[1:].split('||')
                labels[parts[1].split()[0].strip()] = parts[2].strip()
    return labels


def load_split1(h5_path, faa_train, faa_val, faa_test):
    """Load split1 (coarse 5-class) data from h5 + faa files."""
    train_lbl = parse_faa_labels(faa_train)
    val_lbl   = parse_faa_labels(faa_val)
    test_lbl  = parse_faa_labels(faa_test)
    coarse2idx = {c: i for i, c in enumerate(COARSE_CLASSES)}

    X = defaultdict(list); y = defaultdict(list)
    print(f"  Loading {h5_path} ...")
    with h5py.File(h5_path, 'r') as f:
        for pid in f:
            emb = f[pid][:]
            for part, lbl_dict in [('train', train_lbl), ('val', val_lbl), ('test', test_lbl)]:
                if pid in lbl_dict:
                    X[part].append(emb)
                    y[part].append(coarse2idx[lbl_dict[pid]])
                    break

    return (np.array(X['train']), np.array(y['train']),
            np.array(X['val']),   np.array(y['val']),
            np.array(X['test']),  np.array(y['test'])), coarse2idx


def load_split2(h5_path, faa_train, faa_val, faa_test):
    """
    Load all split2 embeddings and labels into memory.
    Returns: {protein_id: (embedding, fine_label, partition)}
    """
    train_lbl = parse_faa_labels(faa_train)
    val_lbl   = parse_faa_labels(faa_val)
    test_lbl  = parse_faa_labels(faa_test)

    data = {}
    print(f"  Loading {h5_path} ...")
    with h5py.File(h5_path, 'r') as f:
        for pid in f:
            emb = f[pid][:]
            if   pid in train_lbl: data[pid] = (emb, train_lbl[pid], 'train')
            elif pid in val_lbl:   data[pid] = (emb, val_lbl[pid],   'val')
            elif pid in test_lbl:  data[pid] = (emb, test_lbl[pid],  'test')
    print(f"  Loaded {len(data)} embeddings")
    return data


# ── Training ──────────────────────────────────────────────────────────────────

def train_model(X_train, y_train, X_val, y_val, n_classes,
                lr=1e-3, batch_size=128, max_epochs=200, patience=10, epsilon=1e-3):
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

    X_tr = torch.tensor(X_train, dtype=torch.float32)
    y_tr = torch.tensor(y_train, dtype=torch.long)
    X_va = torch.tensor(X_val,   dtype=torch.float32).to(device)
    y_va = torch.tensor(y_val,   dtype=torch.long).to(device)

    loader    = DataLoader(TensorDataset(X_tr, y_tr), batch_size=batch_size, shuffle=True,
                           generator=torch.Generator().manual_seed(SEED))
    model     = FNN(n_features=X_train.shape[1], n_classes=n_classes).to(device)
    optimizer = optim.Adam(model.parameters(), lr=lr)
    criterion = nn.CrossEntropyLoss()

    best_val, best_state, no_improve = float('inf'), None, 0

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

        if (epoch + 1) % 10 == 0:
            print(f"    Epoch {epoch+1:3d}  val_loss={val_loss:.4f}")

        if val_loss < best_val - epsilon:
            best_val  = val_loss
            best_state = {k: v.cpu().clone() for k, v in model.state_dict().items()}
            no_improve = 0
        else:
            no_improve += 1
            if no_improve >= patience:
                print(f"    Early stopping at epoch {epoch+1}  (best val_loss={best_val:.4f})")
                break

    model.load_state_dict(best_state)
    return model.to(device), device


def predict_model(model, X, device):
    model.eval()
    with torch.no_grad():
        return model(torch.tensor(X, dtype=torch.float32).to(device)).argmax(1).cpu().numpy()


# ── Hierarchical inference ────────────────────────────────────────────────────

def predict_hierarchical(l1_model, sub_models, X_test, true_fine_lbls, device,
                         coarse2idx, oracle=False):
    """
    Predict fine-grained labels hierarchically.
    If oracle=True, uses the true coarse label for routing (upper bound).
    Returns list of predicted fine-grained label strings.
    """
    idx2coarse = {v: k for k, v in coarse2idx.items()}

    if oracle:
        coarse_preds = [COARSE_MAP[lbl] for lbl in true_fine_lbls]
    else:
        l1_idxs      = predict_model(l1_model, X_test, device)
        coarse_preds = [idx2coarse[i] for i in l1_idxs]

    fine_preds = []
    for x, coarse_label in zip(X_test, coarse_preds):
        subclasses = SUBGROUPS[coarse_label]
        if len(subclasses) == 1:
            fine_preds.append(subclasses[0])
        else:
            sub_model, local_labels = sub_models[coarse_label]
            emb = torch.tensor(x, dtype=torch.float32).unsqueeze(0).to(device)
            sub_model.eval()
            with torch.no_grad():
                sub_idx = sub_model(emb).argmax(1).item()
            fine_preds.append(local_labels[sub_idx])

    return fine_preds


# ── Metrics ───────────────────────────────────────────────────────────────────

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
    rng = np.random.default_rng(SEED)
    samples = defaultdict(list)
    for _ in range(n_iter):
        idx = rng.choice(len(y_true), size=len(y_true), replace=True)
        for k, v in compute_metrics(y_true[idx], y_pred[idx]).items():
            samples[k].append(v)
    return {
        k: {'mean': float(np.mean(v)),
            'lower': float(np.percentile(v, 100 * confidence / 2)),
            'upper': float(np.percentile(v, 100 * (1 - confidence / 2)))}
        for k, v in samples.items()
    }


def per_class_metrics(y_true, y_pred, class_names):
    out = {}
    for i, name in enumerate(class_names):
        mask = y_true == i
        if mask.sum() == 0:
            continue
        out[name] = {
            'accuracy':  round(float(accuracy_score(y_true[mask], y_pred[mask])), 4),
            'precision': round(float(precision_score(y_true == i, y_pred == i, zero_division=0)), 4),
            'recall':    round(float(recall_score(y_true == i, y_pred == i, zero_division=0)), 4),
            'f1_score':  round(float(f1_score(y_true == i, y_pred == i, zero_division=0)), 4),
        }
    return out


def bootstrap_per_class_metrics(y_true, y_pred, class_names, n_iter=30, confidence=0.05):
    rng     = np.random.default_rng(SEED)
    samples = {name: {'accuracy': [], 'precision': [], 'recall': [], 'f1_score': []}
               for name in class_names}
    for _ in range(n_iter):
        idx = rng.choice(len(y_true), size=len(y_true), replace=True)
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
        return {'mean':  float(np.mean(vals)),
                'lower': float(np.percentile(vals, 100 * confidence / 2)),
                'upper': float(np.percentile(vals, 100 * (1 - confidence / 2)))}

    return {name: {metric: ci(vals) for metric, vals in metrics.items() if vals}
            for name, metrics in samples.items()}


# ── Main ──────────────────────────────────────────────────────────────────────

if __name__ == "__main__":

    RESULTS_DIR = "results/training_results"
    os.makedirs(RESULTS_DIR, exist_ok=True)

    # ── Level 1: coarse classifier (trained on split1_30pct) ──────────────────
    print("\n" + "="*60)
    print("  Level 1: Coarse classifier  [split1_30pct, 5 classes]")
    print("="*60)

    (X1_tr, y1_tr, X1_va, y1_va, X1_te, y1_te), coarse2idx = load_split1(
        "embeddings/split1_30pct.h5",
        "dataset_splits/rep_split1_train.faa",
        "dataset_splits/rep_split1_val.faa",
        "dataset_splits/rep_split1_test.faa",
    )
    print(f"  train={len(X1_tr)}  val={len(X1_va)}  test={len(X1_te)}")
    print("  Training ...")
    l1_model, device = train_model(X1_tr, y1_tr, X1_va, y1_va, n_classes=5)
    l1_acc = float(accuracy_score(y1_te, predict_model(l1_model, X1_te, device)))
    print(f"  Level 1 accuracy (split1 test): {l1_acc:.4f}")

    # ── Level 2: sub-classifiers (trained on split2_30pct subsets) ────────────
    print("\n" + "="*60)
    print("  Loading split2_30pct data for Level 2 sub-classifiers")
    print("="*60)

    split2_data = load_split2(
        "embeddings/split2_30pct.h5",
        "dataset_splits/rep_split2_train.faa",
        "dataset_splits/rep_split2_val.faa",
        "dataset_splits/rep_split2_test.faa",
    )

    sub_models = {}  # {group: (model, local_labels)}
    sub_results = {}

    for group, subclasses in SUBGROUPS.items():
        if len(subclasses) == 1:
            print(f"\n  {group}: single class — no sub-classifier needed.")
            continue

        local_labels = subclasses  # already sorted
        loc2idx = {lbl: i for i, lbl in enumerate(local_labels)}

        print(f"\n  Training {group} sub-classifier  ({len(local_labels)} classes) ...")

        X_parts = defaultdict(list); y_parts = defaultdict(list)
        for pid, (emb, fine_lbl, part) in split2_data.items():
            if COARSE_MAP.get(fine_lbl) == group:
                X_parts[part].append(emb)
                y_parts[part].append(loc2idx[fine_lbl])

        X_tr = np.array(X_parts['train']); y_tr = np.array(y_parts['train'])
        X_va = np.array(X_parts['val']);   y_va = np.array(y_parts['val'])
        X_te = np.array(X_parts['test']);  y_te = np.array(y_parts['test'])
        print(f"    train={len(X_tr)}  val={len(X_va)}  test={len(X_te)}")

        model, dev = train_model(X_tr, y_tr, X_va, y_va, n_classes=len(local_labels))
        sub_preds  = predict_model(model, X_te, dev)
        sub_acc    = float(accuracy_score(y_te, sub_preds))
        sub_f1     = float(f1_score(y_te, sub_preds, average='macro', zero_division=0))
        print(f"    {group} sub-classifier: accuracy={sub_acc:.4f}  macro-F1={sub_f1:.4f}")

        sub_models[group]  = (model, local_labels)
        sub_results[group] = {'accuracy': sub_acc, 'macro_f1': sub_f1,
                              'classes': local_labels}

    # ── Hierarchical evaluation on split2_30pct test set ──────────────────────
    print("\n" + "="*60)
    print("  Hierarchical evaluation  [split2_30pct test set, 12 classes]")
    print("="*60)

    fine2idx = {lbl: i for i, lbl in enumerate(FINE_CLASSES)}

    # Build test set
    X_test, y_test_lbls = [], []
    for pid, (emb, fine_lbl, part) in split2_data.items():
        if part == 'test':
            X_test.append(emb)
            y_test_lbls.append(fine_lbl)
    X_test     = np.array(X_test)
    y_test_true = np.array([fine2idx[lbl] for lbl in y_test_lbls])
    print(f"  Test set: {len(X_test)} samples")

    # Hierarchical prediction
    pred_lbls  = predict_hierarchical(l1_model, sub_models, X_test, y_test_lbls,
                                      device, coarse2idx, oracle=False)
    y_test_pred = np.array([fine2idx[lbl] for lbl in pred_lbls])

    # Oracle prediction (upper bound: perfect Level 1 routing)
    oracle_lbls  = predict_hierarchical(l1_model, sub_models, X_test, y_test_lbls,
                                        device, coarse2idx, oracle=True)
    y_oracle_pred = np.array([fine2idx[lbl] for lbl in oracle_lbls])

    # Metrics
    metrics       = compute_metrics(y_test_true, y_test_pred)
    oracle_metrics = compute_metrics(y_test_true, y_oracle_pred)
    boot          = bootstrap_metrics(y_test_true, y_test_pred)
    per_cls       = per_class_metrics(y_test_true, y_test_pred, FINE_CLASSES)
    print("  Running per-class bootstrap (n=30) ...")
    per_cls_boot  = bootstrap_per_class_metrics(y_test_true, y_test_pred, FINE_CLASSES)

    print(f"\n  ── Hierarchical test results ──")
    for k, v in metrics.items():
        print(f"    {k:28s}: {v:.4f}")
    print(f"\n  ── Oracle results (perfect Level 1) ──")
    print(f"    accuracy                    : {oracle_metrics['accuracy']:.4f}")
    print(f"    matthews-corr-coeff         : {oracle_metrics['matthews-corr-coeff']:.4f}")
    print(f"\n  ── Per-class F1 (hierarchical) ──")
    for cls, m in per_cls.items():
        print(f"    {cls:<25s}: f1={m['f1_score']:.3f}")

    # ── Save ──────────────────────────────────────────────────────────────────
    result = {
        'approach':   'hierarchical',
        'split':      'split2_30pct',
        'embedder':   'ProtT5',
        'level1': {
            'training_split': 'split1_30pct',
            'classes':        COARSE_CLASSES,
            'test_accuracy':  l1_acc,
        },
        'level2_sub_classifiers': sub_results,
        'n_test':                 int(len(X_test)),
        'classes':                FINE_CLASSES,
        'test_metrics':                metrics,
        'test_metrics_bootstrap':      boot,
        'oracle_metrics':              oracle_metrics,
        'per_class_metrics':           per_cls,
        'per_class_metrics_bootstrap': per_cls_boot,
    }

    out_path = os.path.join(RESULTS_DIR, "split2_30pct_hierarchical_result.json")
    with open(out_path, 'w') as f:
        json.dump(result, f, indent=2)
    print(f"\n  Saved → {out_path}")
