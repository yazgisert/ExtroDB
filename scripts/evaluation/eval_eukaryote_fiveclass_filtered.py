# 5-class ProtT5 model on the leakage-filtered algae: accuracy, macro-F1 and confusion matrix.

import json

import h5py
import numpy as np
import torch
import torch.nn as nn
from sklearn.metrics import accuracy_score, confusion_matrix, f1_score

ALGAE_H5     = "embeddings/eukaryote_algae.h5"
ALGAE_LABELS = "embeddings/eukaryote_algae_labels.tsv"
BLACKLIST    = "results/eukaryote_leakage/flagged_protein_ids.txt"
MODEL_PATH   = "results/training_results/fiveclass_extremophile_model.pt"
OUT_JSON     = "results/training_results/eukaryote_fiveclass_filtered_result.json"

SEED = 42
N_BOOTSTRAP = 1000
CLASS_NAMES = ["acidophile", "alkaliphile", "halophile", "psychrophile", "thermophile"]


class FNN(nn.Module):
    def __init__(self, n_features=1024, n_classes=5, dropout_rate=0.25):
        super().__init__()
        self.fc1     = nn.Linear(n_features, 32)
        self.dropout = nn.Dropout(dropout_rate)
        self.fc2     = nn.Linear(32, n_classes)

    def forward(self, x):
        return self.fc2(self.dropout(torch.relu(self.fc1(x))))


if __name__ == "__main__":
    with open(BLACKLIST) as f:
        blacklist = {line.strip() for line in f if line.strip()}

    labels = {}
    with open(ALGAE_LABELS) as f:
        next(f)
        for line in f:
            protein_id, _, broad_class = line.rstrip("\n").split("\t")
            if protein_id not in blacklist and broad_class in CLASS_NAMES:
                labels[protein_id] = broad_class

    ids, X = [], []
    with h5py.File(ALGAE_H5, "r") as f:
        for seq_id in f:
            if seq_id in labels:
                ids.append(seq_id)
                X.append(f[seq_id][:])
    X = np.array(X, dtype=np.float32)
    y_true = np.array([CLASS_NAMES.index(labels[i]) for i in ids])
    print(f"{len(ids):,} filtered algae proteins")

    model = FNN()
    model.load_state_dict(torch.load(MODEL_PATH, map_location="cpu"))
    model.eval()
    with torch.no_grad():
        y_pred = model(torch.from_numpy(X)).argmax(dim=1).numpy()

    present = sorted(set(y_true))
    cm = confusion_matrix(y_true, y_pred, labels=range(len(CLASS_NAMES)))
    cm_norm = cm / np.maximum(cm.sum(axis=1, keepdims=True), 1)

    print("\n" + "true / predicted".ljust(16) + "".join(f"{c[:6]:>9}" for c in CLASS_NAMES))
    for i in present:
        print(f"  {CLASS_NAMES[i]:<14}" + "".join(f"{v:>9.3f}" for v in cm_norm[i]))
    accuracy = accuracy_score(y_true, y_pred)
    macro_f1 = f1_score(y_true, y_pred, labels=present, average="macro", zero_division=0)
    print(f"\n  accuracy {accuracy:.3f}   macro-F1 (4 algae classes) {macro_f1:.3f}")

    point = f1_score(y_true, y_pred, labels=present, average=None, zero_division=0)
    rng = np.random.default_rng(SEED)
    boot = []
    for _ in range(N_BOOTSTRAP):
        idx = rng.choice(len(y_true), len(y_true), replace=True)
        boot.append(f1_score(y_true[idx], y_pred[idx], labels=present, average=None, zero_division=0))
    boot = np.array(boot)
    per_class_f1 = {CLASS_NAMES[c]: {"value": float(point[k]),
                                     "ci_lower": float(np.percentile(boot[:, k], 2.5)),
                                     "ci_upper": float(np.percentile(boot[:, k], 97.5))}
                    for k, c in enumerate(present)}

    with open(OUT_JSON, "w") as f:
        json.dump({
            "task":    "5-class model on leakage-filtered eukaryote algae",
            "model":   MODEL_PATH,
            "classes": CLASS_NAMES,
            "eukaryote_eval": {
                "n": len(ids),
                "class_counts": {c: int((y_true == i).sum()) for i, c in enumerate(CLASS_NAMES)},
                "accuracy": float(accuracy),
                "macro_f1": float(macro_f1),
                "per_class_f1": per_class_f1,
                "confusion_matrix": cm.tolist(),
                "confusion_matrix_normalised": cm_norm.tolist(),
            },
        }, f, indent=2)
    print(f"\nResults saved to {OUT_JSON}")
