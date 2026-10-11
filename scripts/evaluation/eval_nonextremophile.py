# Runs the 5-class model (trained on extremophiles only) on a random non-extremophile sample as large as the
# test set. The model has no non-extremophile class, so this reports which class these proteins get and how
# confident the model is.

import json
import os
import random
from collections import defaultdict, Counter

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

EXTREM_H5        = "embeddings/split1_30pct.h5"
EXTREM_TRAIN_FAA = "dataset_splits/rep_split1_train.faa"
EXTREM_VAL_FAA   = "dataset_splits/rep_split1_val.faa"
EXTREM_TEST_FAA  = "dataset_splits/rep_split1_test.faa"

NONEXTREM_H5     = "embeddings/nonextremophile.h5"

RESULTS_DIR = "results/training_results"
OUT_JSON    = os.path.join(RESULTS_DIR, "nonextremophile_eval_result.json")

CLASS_NAMES = ['acidophile', 'alkaliphile', 'halophile', 'psychrophile', 'thermophile']
LABEL2IDX   = {c: i for i, c in enumerate(CLASS_NAMES)}


class FNN(nn.Module):
    def __init__(self, n_features=1024, n_classes=5, dropout_rate=0.25):
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
            if line.startswith('>'):
                parts = line[1:].split('||')
                protein_id = parts[1].split()[0].strip()
                label      = parts[2].strip()
                labels[protein_id] = label
    return labels


def load_extremophile_split(h5_path, faa_train, faa_val, faa_test):
    train_labels = parse_faa_labels(faa_train)
    val_labels   = parse_faa_labels(faa_val)
    test_labels  = parse_faa_labels(faa_test)

    all_labels = sorted(set(train_labels.values()) |
                        set(val_labels.values())   |
                        set(test_labels.values()))
    label2idx  = {lbl: i for i, lbl in enumerate(all_labels)}

    X = defaultdict(list)
    y = defaultdict(list)

    print(f"  Reading extremophile embeddings from {h5_path} ...")
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


def load_nonextremophile_sample(h5_path, n_sample, rng):
    print(f"  Reading non-extremophile embeddings from {h5_path} ...")
    all_embs = []
    with h5py.File(h5_path, 'r') as f:
        keys = list(f.keys())
        print(f"  Total available: {len(keys):,}, sampling {n_sample:,}")
        chosen = rng.choice(len(keys), size=n_sample, replace=False)
        for idx in chosen:
            all_embs.append(f[keys[idx]][:])
    return np.array(all_embs)


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

    best_val_loss = float('inf')
    best_state    = None
    no_improve    = 0
    history       = {'train_loss': [], 'val_loss': []}

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
            print(f"  Epoch {epoch+1:3d}  train={train_loss:.4f}  val={val_loss:.4f}")

        if val_loss < best_val_loss - epsilon:
            best_val_loss = val_loss
            best_state    = {k: v.cpu().clone() for k, v in model.state_dict().items()}
            no_improve    = 0
        else:
            no_improve += 1
            if no_improve >= patience:
                print(f"  Early stopping at epoch {epoch+1}  (best val: {best_val_loss:.4f})")
                break

    model.load_state_dict(best_state)
    return model.to(device), history, device


def predict_with_probs(model, X, device):
    model.eval()
    with torch.no_grad():
        logits = model(torch.tensor(X, dtype=torch.float32).to(device))
        probs  = torch.softmax(logits, dim=1).cpu().numpy()
    preds = probs.argmax(axis=1)
    return preds, probs


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


def main():
    rng = np.random.default_rng(SEED)

    print("\nStep 1/4: loading extremophile data")
    X_train, y_train, X_val, y_val, X_test, y_test, class_names = \
        load_extremophile_split(EXTREM_H5, EXTREM_TRAIN_FAA,
                                EXTREM_VAL_FAA, EXTREM_TEST_FAA)
    n_classes = len(class_names)
    print(f"  train={len(X_train):,}  val={len(X_val):,}  test={len(X_test):,}")
    print(f"  classes: {class_names}")

    print("\nStep 2/4: loading non-extremophile sample")
    X_nonext = load_nonextremophile_sample(NONEXTREM_H5, n_sample=len(X_test), rng=rng)

    print("\nStep 3/4: training the 5-class FNN")
    model, history, device = train(X_train, y_train, X_val, y_val, n_classes=n_classes)

    print("\nStep 4/4: evaluating")

    y_pred_ext, probs_ext = predict_with_probs(model, X_test, device)
    metrics_ext   = compute_metrics(y_test, y_pred_ext)
    per_cls_ext   = per_class_metrics(y_test, y_pred_ext, class_names)
    max_prob_ext  = probs_ext.max(axis=1).tolist()

    print("\n  Extremophile test set:")
    for k, v in metrics_ext.items():
        print(f"    {k:28s}: {v:.4f}")
    print(f"    Mean max-softmax prob       : {np.mean(max_prob_ext):.4f}")

    y_pred_nxt, probs_nxt = predict_with_probs(model, X_nonext, device)
    max_prob_nxt  = probs_nxt.max(axis=1).tolist()
    pred_counts   = Counter(int(p) for p in y_pred_nxt)
    pred_dist     = {class_names[i]: pred_counts.get(i, 0) for i in range(n_classes)}
    pred_dist_pct = {k: round(v / len(y_pred_nxt) * 100, 2) for k, v in pred_dist.items()}

    print("\n  Non-extremophile proteins (n={:,}):".format(len(X_nonext)))
    print(f"    Mean max-softmax prob       : {np.mean(max_prob_nxt):.4f}")
    print(f"    Predicted class distribution:")
    for cls, cnt in pred_dist.items():
        print(f"      {cls:20s}: {cnt:6,}  ({pred_dist_pct[cls]:.1f}%)")

    print("\n  Confidence comparison:")
    print(f"    Extremophile  mean max-prob : {np.mean(max_prob_ext):.4f}")
    print(f"    Non-extremoph mean max-prob : {np.mean(max_prob_nxt):.4f}")
    print(f"    Gap                         : {np.mean(max_prob_ext) - np.mean(max_prob_nxt):.4f}")

    os.makedirs(RESULTS_DIR, exist_ok=True)
    result = {
        "experiment":          "nonextremophile_eval",
        "split":               "split1_30pct",
        "class_names":         class_names,
        "n_train":             int(len(X_train)),
        "n_val":               int(len(X_val)),
        "n_test_extremophile": int(len(X_test)),
        "n_test_nonextreme":   int(len(X_nonext)),
        "note": (
            "Non-extremophile proteins sampled from nonextremophile.h5 "
            "(UniProt, deduped against extremophile set at 90% identity). "
            "The model was never trained on non-extremophile proteins and "
            "cannot abstain — it always predicts one of the 5 extremophile classes."
        ),

        "extremophile_test_metrics": metrics_ext,
        "extremophile_per_class":    per_cls_ext,
        "extremophile_max_prob":     [round(p, 4) for p in max_prob_ext],
        "extremophile_mean_max_prob": round(float(np.mean(max_prob_ext)), 4),

        "nonextreme_predicted_class_counts":   pred_dist,
        "nonextreme_predicted_class_pct":      pred_dist_pct,
        "nonextreme_max_prob":        [round(p, 4) for p in max_prob_nxt],
        "nonextreme_mean_max_prob":   round(float(np.mean(max_prob_nxt)), 4),
        "nonextreme_predicted_labels": [class_names[int(p)] for p in y_pred_nxt],

        # Per-class softmax probabilities (for detailed plotting)
        "extremophile_all_probs":  [[round(float(v), 4) for v in row] for row in probs_ext],
        "nonextreme_all_probs":    [[round(float(v), 4) for v in row] for row in probs_nxt],
    }

    with open(OUT_JSON, 'w') as f:
        json.dump(result, f, indent=2)
    print(f"\nResults saved to {OUT_JSON}")


if __name__ == "__main__":
    main()
