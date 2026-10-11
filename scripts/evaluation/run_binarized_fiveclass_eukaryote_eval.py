# Binarized 5-class ProtT5 model on the leakage-filtered algae: recall overall and per class, printed next to the
# dedicated binary model.

import json
import os
import random
from collections import defaultdict

import h5py
import numpy as np
import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import DataLoader, TensorDataset

SEED = 42
random.seed(SEED)
np.random.seed(SEED)
torch.manual_seed(SEED)

EXTREM_H5        = "embeddings/split1_30pct.h5"
EXTREM_TRAIN_FAA = "dataset_splits/rep_split1_train.faa"
EXTREM_VAL_FAA   = "dataset_splits/rep_split1_val.faa"
EXTREM_TEST_FAA  = "dataset_splits/rep_split1_test.faa"
ALGAE_H5         = "embeddings/eukaryote_algae.h5"
ALGAE_LABELS_TSV = "embeddings/eukaryote_algae_labels.tsv"
BLACKLIST_PATH   = "results/eukaryote_leakage/flagged_protein_ids.txt"
MODEL_PATH       = "results/training_results/fiveclass_extremophile_model.pt"
BINARY_JSON      = "results/training_results/binary_eukaryote_eval_clustered_result.json"
OUT_JSON         = "results/training_results/binarized_fiveclass_eukaryote_eval_result.json"

# Sorted alphabetically — must match the order used during training
CLASS_NAMES = ["acidophile", "alkaliphile", "halophile", "psychrophile", "thermophile"]
LABEL2IDX   = {c: i for i, c in enumerate(CLASS_NAMES)}

BINARIZATION_THRESHOLD = 0.5


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
            if line.startswith(">"):
                parts = line[1:].split("||")
                protein_id = parts[1].split()[0].strip()
                labels[protein_id] = parts[2].strip()
    return labels


def load_prokaryote_splits(h5_path, faa_train, faa_val, faa_test):
    train_labels = parse_faa_labels(faa_train)
    val_labels   = parse_faa_labels(faa_val)
    test_labels  = parse_faa_labels(faa_test)

    X = defaultdict(list)
    y = defaultdict(list)

    with h5py.File(h5_path, "r") as f:
        for seq_id in f:
            emb = f[seq_id][:]
            for split, label_map in [("train", train_labels), ("val", val_labels), ("test", test_labels)]:
                if seq_id in label_map:
                    lbl = label_map[seq_id]
                    if lbl in LABEL2IDX:
                        X[split].append(emb)
                        y[split].append(LABEL2IDX[lbl])
                    break

    return (np.array(X["train"]), np.array(y["train"]),
            np.array(X["val"]),   np.array(y["val"]),
            np.array(X["test"]),  np.array(y["test"]))


def train_model(X_train, y_train, X_val, y_val,
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

    model     = FNN(n_features=X_train.shape[1], n_classes=5).to(device)
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
    softmax   = nn.Softmax(dim=1)
    all_probs = []

    with torch.no_grad():
        for i in range(0, len(X), batch_size):
            batch = torch.tensor(X[i:i + batch_size], dtype=torch.float32).to(device)
            all_probs.append(softmax(model(batch)).cpu().numpy())

    return np.concatenate(all_probs)


if __name__ == "__main__":

    print("Loading prokaryote training data ...")
    X_train, y_train, X_val, y_val, X_test, y_test = load_prokaryote_splits(
        EXTREM_H5, EXTREM_TRAIN_FAA, EXTREM_VAL_FAA, EXTREM_TEST_FAA
    )
    print(f"  train: {len(X_train):,}  val: {len(X_val):,}  test: {len(X_test):,}")

    print("\nTraining 5-class model ...")
    model, device = train_model(X_train, y_train, X_val, y_val)

    os.makedirs(os.path.dirname(MODEL_PATH), exist_ok=True)
    torch.save(model.state_dict(), MODEL_PATH)
    print(f"  Saved to {MODEL_PATH}")

    print("\nLoading blacklist ...")
    with open(BLACKLIST_PATH) as f:
        blacklist = set(line.strip() for line in f if line.strip())
    print(f"  {len(blacklist):,} proteins will be skipped")

    print("\nLoading filtered eukaryote labels ...")
    protein_to_class = {}
    with open(ALGAE_LABELS_TSV) as f:
        next(f)
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) == 3:
                pid, _, broad_class = parts
                if pid not in blacklist and broad_class in LABEL2IDX:
                    protein_to_class[pid] = broad_class

    print("\nLoading filtered eukaryote embeddings ...")
    algae_ids    = []
    algae_embeds = []
    with h5py.File(ALGAE_H5, "r") as f:
        for seq_id in f:
            if seq_id in protein_to_class:
                algae_ids.append(seq_id)
                algae_embeds.append(f[seq_id][:])

    X_algae = np.array(algae_embeds)
    print(f"  {len(X_algae):,} proteins after filtering")

    print("\nRunning 5-class inference ...")
    probs = get_probs(model, X_algae, device)  # shape: (n_proteins, 5)

    # Binarize: max probability across all 5 classes must reach the threshold
    max_probs      = probs.max(axis=1)
    binary_preds   = (max_probs >= BINARIZATION_THRESHOLD).astype(int)

    overall_recall = float(binary_preds.mean())
    print(f"\n  Overall recall (binarized, threshold={BINARIZATION_THRESHOLD}): {overall_recall:.4f}")

    class_buckets = defaultdict(list)
    for i, pid in enumerate(algae_ids):
        class_buckets[protein_to_class[pid]].append(binary_preds[i])

    # Dedicated binary model (clustered split) on the same filtered algae, for the printed comparison
    binary_recalls = {}
    if os.path.exists(BINARY_JSON):
        with open(BINARY_JSON) as f:
            binary_recalls = json.load(f).get("by_class", {})

    print(f"\n  Per-class recall (binarized 5-class vs. dedicated binary model):")
    print(f"  {'Class':<15} {'n':>8} {'Binarized 5-class':>20} {'Dedicated binary':>18}")
    print("  " + "-" * 65)

    class_summary = {}
    for cls in sorted(class_buckets):
        preds  = np.array(class_buckets[cls])
        recall = float(preds.mean())
        n      = len(preds)
        bin_r  = binary_recalls.get(cls, {}).get("recall")
        bin_str = f"{bin_r:.3f}" if bin_r is not None else "n/a"
        print(f"  {cls:<15} {n:>8,} {recall:>20.3f} {bin_str:>18}")
        class_summary[cls] = {"n_proteins": n, "recall": recall}

    result = {
        "task":                "binarized_5class_on_filtered_eukaryote_algae",
        "model":               "ProtT5 FNN 5-class (split1_30pct), binarized by max-probability threshold",
        "binarization_threshold": BINARIZATION_THRESHOLD,
        "filtering":           f"Removed proteins >= 30% identity to prokaryote train+val (blacklist: {len(blacklist):,} proteins)",
        "n_algae_proteins":    int(len(X_algae)),
        "overall_recall":      overall_recall,
        "by_class":            class_summary,
    }

    with open(OUT_JSON, "w") as f:
        json.dump(result, f, indent=2)
    print(f"\nResults saved to {OUT_JSON}")
