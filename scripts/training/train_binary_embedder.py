# Extremophile vs. non-extremophile FNN on the clustered binary split, for one embedder (--embedder prott5 / esm2_650m
# / esm2_8m / onehot). Saves the model and test metrics with 1000 bootstrap CIs.

import argparse
import json
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

SPLIT_TSV = "dataset_splits/clustered_binary/binary_split.tsv"
EXTREM_FAA = [
    "dataset_splits/rep_split1_train.faa",
    "dataset_splits/rep_split1_val.faa",
    "dataset_splits/rep_split1_test.faa",
]
NONEXTREM_FAA = "nonextremophile_data/nonextremophile_filtered.faa"
RESULTS_DIR   = "results/training_results"

EMBEDDINGS = {
    "prott5":    ("embeddings/split1_30pct.h5",        "embeddings/nonextremophile.h5"),
    "esm2_650m": ("embeddings/split1_30pct_esm2.h5",   "embeddings/nonextremophile_esm2.h5"),
    "esm2_8m":   ("embeddings/split1_30pct_esm2_8m.h5", "embeddings/nonextremophile_esm2_8m.h5"),
    "onehot":    None,
}

AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"
N_BOOTSTRAP = 1000


class FNN(nn.Module):
    def __init__(self, n_features, n_classes=2, dropout_rate=0.25):
        super().__init__()
        self.fc1     = nn.Linear(n_features, 32)
        self.dropout = nn.Dropout(dropout_rate)
        self.fc2     = nn.Linear(32, n_classes)

    def forward(self, x):
        return self.fc2(self.dropout(torch.relu(self.fc1(x))))


def read_faa(filepath, id_from_header):
    seqs = {}
    current_id, current_seq = None, []
    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if current_id is not None:
                    seqs[current_id] = "".join(current_seq)
                current_id, current_seq = id_from_header(line[1:]), []
            else:
                current_seq.append(line)
    if current_id is not None:
        seqs[current_id] = "".join(current_seq)
    return seqs


def split_id(header):
    return header.split("||")[1].split()[0].strip()


def amino_acid_composition(seq):
    counts = np.array([seq.upper().count(aa) for aa in AMINO_ACIDS], dtype=np.float32)
    return counts / max(counts.sum(), 1)


def read_split():
    extrem, nonextrem = defaultdict(list), defaultdict(list)
    with open(SPLIT_TSV) as f:
        next(f)
        for line in f:
            protein_id, split, label = line.rstrip("\n").split("\t")
            (nonextrem if label == "non_extremophile" else extrem)[split].append(protein_id)
    return extrem, nonextrem


def load_features(embedder, extrem, nonextrem):
    if embedder == "onehot":
        seqs = {}
        for path in EXTREM_FAA:
            seqs.update(read_faa(path, split_id))
        seqs.update(read_faa(NONEXTREM_FAA, lambda h: h.split()[0]))
        lookup = lambda ids, _: [amino_acid_composition(seqs[i]) for i in ids]
        sources = (None, None)
    else:
        def lookup(ids, path):
            with h5py.File(path, "r") as f:
                missing = [i for i in ids if i not in f]
                if missing:
                    raise SystemExit(f"{len(missing):,} proteins missing from {path}, e.g. {missing[:3]}. "
                                     f"Finish the embeddings first.")
                return [f[i][:] for i in ids]
        sources = EMBEDDINGS[embedder]

    rng = np.random.default_rng(SEED)
    data = {}
    for split in ("train", "val", "test"):
        X = np.array(lookup(extrem[split], sources[0]) + lookup(nonextrem[split], sources[1]), dtype=np.float32)
        y = np.array([1] * len(extrem[split]) + [0] * len(nonextrem[split]))
        idx = rng.permutation(len(X))
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

    model     = FNN(n_features=X_train.shape[1]).to(device)
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


def binary_metrics(y_true, y_pred):
    return {
        "accuracy":            accuracy_score(y_true, y_pred),
        "balanced_accuracy":   balanced_accuracy_score(y_true, y_pred),
        "precision":           precision_score(y_true, y_pred, zero_division=0),
        "recall":              recall_score(y_true, y_pred, zero_division=0),
        "specificity":         recall_score(y_true, y_pred, pos_label=0, zero_division=0),
        "f1_extremophile":     f1_score(y_true, y_pred, zero_division=0),
        "f1_non_extremophile": f1_score(y_true, y_pred, pos_label=0, zero_division=0),
        "macro_f1":            f1_score(y_true, y_pred, average="macro", zero_division=0),
        "mcc":                 matthews_corrcoef(y_true, y_pred),
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
    parser = argparse.ArgumentParser()
    parser.add_argument("--embedder", choices=EMBEDDINGS, required=True)
    embedder = parser.parse_args().embedder

    print("Reading the clustered binary split ...")
    extrem, nonextrem = read_split()
    for split in ("train", "val", "test"):
        print(f"  {split:<5} extremophile: {len(extrem[split]):,}  non-extremophile: {len(nonextrem[split]):,}")

    print(f"\nLoading {embedder} features ...")
    data = load_features(embedder, extrem, nonextrem)
    X_train, y_train = data["train"]
    X_val,   y_val   = data["val"]
    X_test,  y_test  = data["test"]
    print(f"  {X_train.shape[1]} features per protein")

    print("\nTraining ...")
    model, device = train_model(X_train, y_train, X_val, y_val)
    model_path = f"{RESULTS_DIR}/binary_{embedder}_clustered_model.pt"
    torch.save(model.state_dict(), model_path)
    print(f"  Saved to {model_path}")

    model.eval()
    with torch.no_grad():
        y_pred = model(torch.from_numpy(X_test).to(device)).argmax(dim=1).cpu().numpy()
    metrics = with_ci(y_test, y_pred)

    print(f"\nTest set ({len(y_test):,} proteins)")
    for k, v in metrics.items():
        print(f"  {k:<20} {v['value']:.3f}  [{v['ci_lower']:.3f}, {v['ci_upper']:.3f}]")

    result = {
        "task":       "binary_extremophile_vs_nonextremophile",
        "embedder":   embedder,
        "split":      SPLIT_TSV + " (13,031 non-extremophile cluster representatives + as many extremophiles, equal per class)",
        "n_features": int(X_train.shape[1]),
        "n_train":    int(len(y_train)),
        "n_val":      int(len(y_val)),
        "n_test":     int(len(y_test)),
        "test_metrics": metrics,
    }
    out_json = f"{RESULTS_DIR}/binary_{embedder}_clustered_result.json"
    with open(out_json, "w") as f:
        json.dump(result, f, indent=2)
    print(f"\nResults saved to {out_json}")
