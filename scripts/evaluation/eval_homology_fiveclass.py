# Homology-based inference on the 5-class test set (best hit among extremophile training proteins), compared with
# ProtT5 overall and on proteins with and without a hit.

import json
from collections import defaultdict

import h5py
import numpy as np
import torch
import torch.nn as nn
from sklearn.metrics import accuracy_score, f1_score

SEED = 42

SEARCH_DIR = "results/test_vs_train_search"
EXTREM_H5  = "embeddings/split1_30pct.h5"
MODEL_PATH = "results/training_results/fiveclass_extremophile_model.pt"
OUT_JSON   = "results/training_results/homology_fiveclass_result.json"

CLASS_NAMES = ["acidophile", "alkaliphile", "halophile", "psychrophile", "thermophile"]
NO_HIT      = len(CLASS_NAMES)
N_BOOTSTRAP = 1000


class FNN(nn.Module):
    def __init__(self, n_features=1024, n_classes=5, dropout_rate=0.25):
        super().__init__()
        self.fc1     = nn.Linear(n_features, 32)
        self.dropout = nn.Dropout(dropout_rate)
        self.fc2     = nn.Linear(32, n_classes)

    def forward(self, x):
        return self.fc2(self.dropout(torch.relu(self.fc1(x))))


def metrics(y_true, y_pred):
    return {
        "accuracy": accuracy_score(y_true, y_pred),
        "macro_f1": f1_score(y_true, y_pred, labels=range(len(CLASS_NAMES)), average="macro", zero_division=0),
    }


def with_ci(y_true, y_pred):
    point = metrics(y_true, y_pred)
    rng   = np.random.default_rng(SEED)
    boot  = defaultdict(list)
    for _ in range(N_BOOTSTRAP):
        idx = rng.choice(len(y_true), len(y_true), replace=True)
        for k, v in metrics(y_true[idx], y_pred[idx]).items():
            boot[k].append(v)
    return {k: {"value": float(v),
                "ci_lower": float(np.percentile(boot[k], 2.5)),
                "ci_upper": float(np.percentile(boot[k], 97.5))}
            for k, v in point.items()}


def accuracy_ci(correct):
    # 95% bootstrap CI of the share of correct predictions in a subset
    rng  = np.random.default_rng(SEED)
    boot = [correct[rng.choice(len(correct), len(correct), replace=True)].mean() for _ in range(N_BOOTSTRAP)]
    return [float(np.percentile(boot, 2.5)), float(np.percentile(boot, 97.5))]

def per_class_recall(y_true, y_pred):
    rng = np.random.default_rng(SEED)
    result = {}
    for i, c in enumerate(CLASS_NAMES):
        correct = (y_pred[y_true == i] == i).astype(float)
        boot = [correct[rng.choice(len(correct), len(correct), replace=True)].mean() for _ in range(N_BOOTSTRAP)]
        result[c] = {"value": float(correct.mean()),
                     "ci_lower": float(np.percentile(boot, 2.5)),
                     "ci_upper": float(np.percentile(boot, 97.5))}
    return result


if __name__ == "__main__":
    label, which = {}, {}
    with open(f"{SEARCH_DIR}/labels.tsv") as f:
        next(f)
        for line in f:
            protein_id, split, cls = line.rstrip("\n").split("\t")
            label[protein_id], which[protein_id] = cls, split

    test_ids = [i for i, s in which.items() if s == "test" and label[i] in CLASS_NAMES]

    # Best hit (highest bit score) among extremophile training proteins only,
    # since the 5-class model never saw non-extremophiles either
    best = {}
    with open(f"{SEARCH_DIR}/hits.tsv") as f:
        for line in f:
            query, target, _, _, bits = line.split("\t")[:5]
            if label.get(query) in CLASS_NAMES and label.get(target) in CLASS_NAMES:
                bits = float(bits)
                if query not in best or bits > best[query][0]:
                    best[query] = (bits, label[target])

    y_true = np.array([CLASS_NAMES.index(label[i]) for i in test_ids])
    y_hbi  = np.array([CLASS_NAMES.index(best[i][1]) if i in best else NO_HIT for i in test_ids])
    has_hit = y_hbi != NO_HIT

    with h5py.File(EXTREM_H5, "r") as f:
        X = np.array([f[i][:] for i in test_ids], dtype=np.float32)
    model = FNN()
    model.load_state_dict(torch.load(MODEL_PATH, map_location="cpu"))
    model.eval()
    with torch.no_grad():
        y_prott5 = model(torch.from_numpy(X)).argmax(dim=1).numpy()

    result = {
        "task":     "homology-based inference vs. ProtT5 on the split1_30pct 5-class test set",
        "method":   "label of the best MMseqs2 hit (bit score) among extremophile train + val proteins; no hit counts as wrong",
        "n_test":   len(test_ids),
        "n_with_hit": int(has_hit.sum()),
        "homology": {**with_ci(y_true, y_hbi), "per_class_recall": per_class_recall(y_true, y_hbi)},
        "prott5":   {**with_ci(y_true, y_prott5), "per_class_recall": per_class_recall(y_true, y_prott5)},
        "subsets": {
            "with_hit":    {"n": int(has_hit.sum()),
                            "homology_accuracy":    float((y_hbi[has_hit] == y_true[has_hit]).mean()),
                            "homology_accuracy_ci": accuracy_ci(y_hbi[has_hit] == y_true[has_hit]),
                            "prott5_accuracy":      float((y_prott5[has_hit] == y_true[has_hit]).mean()),
                            "prott5_accuracy_ci":   accuracy_ci(y_prott5[has_hit] == y_true[has_hit])},
            "without_hit": {"n": int((~has_hit).sum()),
                            "prott5_accuracy":      float((y_prott5[~has_hit] == y_true[~has_hit]).mean()),
                            "prott5_accuracy_ci":   accuracy_ci(y_prott5[~has_hit] == y_true[~has_hit])},
        },
    }

    print(f"{len(test_ids):,} extremophile test proteins, {has_hit.sum():,} ({has_hit.mean():.1%}) with a hit to an extremophile training protein\n")
    print(f"{'':22}{'Homology':>10}{'ProtT5':>10}")
    for k in ("accuracy", "macro_f1"):
        print(f"  {k:<20}{result['homology'][k]['value']:>10.3f}{result['prott5'][k]['value']:>10.3f}")
    for c in CLASS_NAMES:
        print(f"  recall {c:<13}{result['homology']['per_class_recall'][c]['value']:>10.3f}{result['prott5']['per_class_recall'][c]['value']:>10.3f}")
    s = result["subsets"]
    print(f"\n  proteins with a hit ({s['with_hit']['n']:,}):    homology {s['with_hit']['homology_accuracy']:.3f}   ProtT5 {s['with_hit']['prott5_accuracy']:.3f}")
    print(f"  proteins without a hit ({s['without_hit']['n']:,}): homology   —     ProtT5 {s['without_hit']['prott5_accuracy']:.3f}")

    with open(OUT_JSON, "w") as f:
        json.dump(result, f, indent=2)
    print(f"\nResults saved to {OUT_JSON}")
