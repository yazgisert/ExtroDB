# Archaea check for the 5-class model (split1 30% test set), ProtT5 or amino acid composition (--embedder): accuracy
# per class for archaeal vs. bacterial proteins (= recall, since each group has one true class), accuracy and macro-F1
# per domain (1000 bootstrap CIs), and halophile / thermophile accuracy per subtype and domain.

import argparse
import csv
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path

import h5py
import numpy as np
import torch
from sklearn.metrics import accuracy_score, f1_score

from eval_eukaryote_fiveclass_filtered import CLASS_NAMES, FNN, MODEL_PATH

TEST_FAA   = "dataset_splits/rep_split1_test.faa"
EMBEDDINGS = "embeddings/split1_30pct.h5"
DOMAIN_TSV = "curation/source_tables/proteome_domain_ncbi.tsv"
METADATA   = "dataset_splits/split_input_metadata.tsv"
OUT_JSON   = {"prott5": "results/training_results/fiveclass_domain_bias_result.json",
              "onehot": "results/training_results/fiveclass_domain_bias_onehot_result.json"}
ONEHOT_MODEL  = "results/training_results/fiveclass_onehot_model.pt"
ONEHOT_RESULT = "results/training_results/split1_30pct_onehot_result.json"

SEED = 42
N_BOOTSTRAP = 1000


def with_ci(y_true, y_pred):
    def metrics(t, p):
        return {"accuracy": accuracy_score(t, p), "macro_f1": f1_score(t, p, average="macro")}
    point = metrics(y_true, y_pred)
    rng = np.random.default_rng(SEED)
    boot = defaultdict(list)
    for _ in range(N_BOOTSTRAP):
        idx = rng.choice(len(y_true), len(y_true), replace=True)
        for k, v in metrics(y_true[idx], y_pred[idx]).items():
            boot[k].append(v)
    return {k: {"value": float(v), "ci_lower": float(np.percentile(boot[k], 2.5)),
                "ci_upper": float(np.percentile(boot[k], 97.5))} for k, v in point.items()}


def onehot_model():
    # The 5-class composition model was not saved when it was trained, so it is retrained here with the original
    # code and seed; the macro-F1 must match the reported result before it is used.
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "training"))
    import train_onehot
    if not Path(ONEHOT_MODEL).exists():
        X_tr, y_tr, X_va, y_va, X_te, y_te, classes = train_onehot.load_split(
            "dataset_splits/rep_split1_train.faa", "dataset_splits/rep_split1_val.faa", TEST_FAA)
        assert classes == CLASS_NAMES
        model, _, device = train_onehot.train(X_tr, y_tr, X_va, y_va, n_classes=5)
        macro_f1 = f1_score(y_te, train_onehot.predict(model, X_te, device), average="macro")
        reported = json.load(open(ONEHOT_RESULT))["test_metrics"]["macro-f1_score"]
        assert abs(macro_f1 - reported) < 1e-9, f"retrained macro-F1 {macro_f1:.4f} differs from reported {reported:.4f}"
        torch.save(model.cpu().state_dict(), ONEHOT_MODEL)
        print(f"Retrained composition model reproduces macro-F1 {macro_f1:.4f}; saved to {ONEHOT_MODEL}")
    model = train_onehot.FNN(n_features=20, n_classes=5)
    model.load_state_dict(torch.load(ONEHOT_MODEL, map_location="cpu"))
    return model, train_onehot.aa_composition


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--embedder", choices=["prott5", "onehot"], default="prott5")
    embedder = parser.parse_args().embedder
    domain_of = {r["proteome_id"]: r["domain"] for r in csv.DictReader(open(DOMAIN_TSV), delimiter="\t")}
    subtype_of = {r["proteome_id"]: r["split2_label"] for r in csv.DictReader(open(METADATA), delimiter="\t")}

    ids, labels, domains, subtypes, seqs = [], [], [], [], {}
    with open(TEST_FAA) as f, h5py.File(EMBEDDINGS, "r") as h5:
        for line in f:
            if line.startswith(">"):
                parts = line[1:].strip().split("||")
                protein_id = parts[1].split()[0]
                seqs[protein_id] = []
                # ProtT5 skips the few test proteins without an embedding; composition uses all of them
                if embedder == "onehot" or protein_id in h5:
                    ids.append(protein_id)
                    labels.append(CLASS_NAMES.index(parts[2]))
                    domains.append(domain_of[parts[0]])
                    subtypes.append(subtype_of[parts[0]])
            else:
                seqs[protein_id].append(line.strip())
        if embedder == "prott5":
            X = np.array([h5[i][:] for i in ids], dtype=np.float32)
    y_true, domains, subtypes = np.array(labels), np.array(domains), np.array(subtypes)

    if embedder == "prott5":
        model = FNN()
        model.load_state_dict(torch.load(MODEL_PATH, map_location="cpu"))
    else:
        model, aa_composition = onehot_model()
        X = np.array([aa_composition("".join(seqs[i])) for i in ids], dtype=np.float32)
    model.eval()
    with torch.no_grad():
        y_pred = model(torch.from_numpy(X)).argmax(dim=1).numpy()

    bacterial = domains == "bacteria"
    result = {
        "n_test": len(ids),
        "n_archaeal": int((~bacterial).sum()),
        "all": with_ci(y_true, y_pred),
        "bacteria_only": with_ci(y_true[bacterial], y_pred[bacterial]),
        "archaea_only": with_ci(y_true[~bacterial], y_pred[~bacterial]),
        "per_class": {},
        "predicted_as_by_domain": {},
    }
    print(f"{len(ids):,} test proteins, {(~bacterial).sum():,} archaeal")
    for subset in ("all", "bacteria_only", "archaea_only"):
        s = result[subset]
        print(f"  {subset:<14} accuracy {s['accuracy']['value']:.3f}  macro-F1 {s['macro_f1']['value']:.3f} "
              f"[{s['macro_f1']['ci_lower']:.3f}, {s['macro_f1']['ci_upper']:.3f}]")

    print(f"\n{'class':<13} {'archaea n':>9} {'acc.':>7} {'bacteria n':>10} {'acc.':>7} {'F1 all':>7} {'F1 bact.':>8}")
    f1_all = f1_score(y_true, y_pred, average=None, labels=range(5))
    f1_bac = f1_score(y_true[bacterial], y_pred[bacterial], average=None, labels=range(5))
    for c, name in enumerate(CLASS_NAMES):
        entry = {"f1_all": float(f1_all[c]), "f1_bacteria_only": float(f1_bac[c])}
        for domain in ("archaea", "bacteria"):
            mask = (y_true == c) & (domains == domain)
            entry[domain] = {"n": int(mask.sum()), "accuracy": float((y_pred[mask] == c).mean()) if mask.any() else None}
        result["per_class"][name] = entry
        a, b = entry["archaea"], entry["bacteria"]
        a_r = f"{a['accuracy']:.3f}" if a["accuracy"] is not None else "–"
        print(f"{name:<13} {a['n']:>9} {a_r:>7} {b['n']:>10} {b['accuracy']:>7.3f} {f1_all[c]:>7.3f} {f1_bac[c]:>8.3f}")

    for domain in ("archaea", "bacteria"):
        counts = Counter(y_pred[domains == domain])
        total = sum(counts.values())
        result["predicted_as_by_domain"][domain] = {CLASS_NAMES[c]: counts[c] / total for c in range(5)}
    print("\nShare of predictions per domain:")
    for domain, shares in result["predicted_as_by_domain"].items():
        print(f"  {domain:<9} " + "  ".join(f"{k} {v:.2f}" for k, v in shares.items()))

    result["accuracy_by_subtype_and_domain"] = {}
    print("\nAccuracy by subtype and domain:")
    for name in ("halophile", "thermophile"):
        c = CLASS_NAMES.index(name)
        for subtype in sorted(set(subtypes[y_true == c])):
            for domain in ("archaea", "bacteria"):
                mask = (y_true == c) & (subtypes == subtype) & (domains == domain)
                if mask.any():
                    acc = float((y_pred[mask] == c).mean())
                    result["accuracy_by_subtype_and_domain"][f"{subtype}|{domain}"] = {"n": int(mask.sum()), "accuracy": acc}
                    print(f"  {subtype:<20} {domain:<9} n={mask.sum():>5}  accuracy {acc:.3f}")

    with open(OUT_JSON[embedder], "w") as f:
        json.dump(result, f, indent=2)
    print(f"\nResults saved to {OUT_JSON[embedder]}")
