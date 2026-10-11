# Binary models of all embedders on the leakage-filtered algae: extremophile recall overall, per class and per species
# with bootstrap CIs; ProtT5 also on the unfiltered algae.

import csv
import json
import os

import h5py
import numpy as np
import torch
import torch.nn as nn

SEED = 42

ALGAE_LABELS = "embeddings/eukaryote_algae_labels.tsv"
ALGAE_FAA    = "results/homology_inference/filtered_eukaryote_query.faa"
BLACKLIST    = "results/eukaryote_leakage/flagged_protein_ids.txt"
RESULTS      = "results/training_results"
OUT_JSON     = f"{RESULTS}/eukaryote_embedder_clustered_result.json"
MASTER_INDEX = "searchable_resource/extremophile_master_index.tsv"

ALGAE_EMBEDDINGS = {
    "onehot":    None,
    "esm2_8m":   "embeddings/eukaryote_algae_filtered_esm2_8m.h5",
    "esm2_650m": "embeddings/eukaryote_algae_filtered_esm2.h5",
    "prott5":    "embeddings/eukaryote_algae.h5",
}
CLASSES     = ["acidophile", "halophile", "psychrophile", "thermophile"]
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


def read_faa(filepath):
    seqs = {}
    current_id, current_seq = None, []
    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if current_id is not None:
                    seqs[current_id] = "".join(current_seq)
                current_id, current_seq = line[1:].split()[0], []
            else:
                current_seq.append(line)
    if current_id is not None:
        seqs[current_id] = "".join(current_seq)
    return seqs


def amino_acid_composition(seq):
    counts = np.array([seq.upper().count(aa) for aa in AMINO_ACIDS], dtype=np.float32)
    return counts / max(counts.sum(), 1)


def algae_proteins(filtered=True):
    blacklist = set()
    if filtered:
        with open(BLACKLIST) as f:
            blacklist = {line.strip() for line in f if line.strip()}
    with h5py.File(ALGAE_EMBEDDINGS["prott5"], "r") as f:
        embedded = set(f.keys())

    labels, proteome = {}, {}
    with open(ALGAE_LABELS) as f:
        next(f)
        for line in f:
            protein_id, proteome_id, broad_class = line.rstrip("\n").split("\t")
            if protein_id in embedded and protein_id not in blacklist and broad_class in CLASSES:
                labels[protein_id] = broad_class
                proteome[protein_id] = proteome_id
    return sorted(labels), labels, proteome


def predict(embedder, ids, model_path):
    X = features(embedder, ids)
    model = FNN(n_features=X.shape[1])
    model.load_state_dict(torch.load(model_path, map_location="cpu"))
    model.eval()
    with torch.no_grad():
        return torch.softmax(model(torch.from_numpy(X)), dim=1)[:, 1].numpy()


def features(embedder, ids):
    if embedder == "onehot":
        seqs = read_faa(ALGAE_FAA)
        return np.array([amino_acid_composition(seqs[i]) for i in ids], dtype=np.float32)
    with h5py.File(ALGAE_EMBEDDINGS[embedder], "r") as f:
        return np.array([f[i][:] for i in ids], dtype=np.float32)


def recall_with_ci(predictions):
    rng  = np.random.default_rng(SEED)
    boot = [predictions[rng.choice(len(predictions), len(predictions), replace=True)].mean()
            for _ in range(N_BOOTSTRAP)]
    return {"n_proteins": int(len(predictions)),
            "recall":     float(predictions.mean()),
            "ci_lower":   float(np.percentile(boot, 2.5)),
            "ci_upper":   float(np.percentile(boot, 97.5))}


if __name__ == "__main__":
    ids, labels, proteome = algae_proteins()
    classes = np.array([labels[i] for i in ids])
    proteomes = np.array([proteome[i] for i in ids])
    with open(MASTER_INDEX) as f:
        species = {row["proteome_id"]: row["taxon_name"] for row in csv.DictReader(f, delimiter="\t")}
    print(f"{len(ids):,} filtered algae proteins\n")

    def summarise(prob, classes, proteomes):
        predicted = (prob >= 0.5).astype(int)
        return {
            "overall":    {**recall_with_ci(predicted), "mean_prob": float(prob.mean())},
            "by_class":   {c: {**recall_with_ci(predicted[classes == c]), "mean_prob": float(prob[classes == c].mean())}
                           for c in CLASSES},
            "by_species": {pid: {"species": species.get(pid, pid), "broad_class": str(classes[proteomes == pid][0]),
                                 "n_proteins": int((proteomes == pid).sum()),
                                 "recall": float(predicted[proteomes == pid].mean()),
                                 "mean_prob": float(prob[proteomes == pid].mean())}
                           for pid in sorted(set(proteomes))},
        }

    results = {}
    for embedder, path in ALGAE_EMBEDDINGS.items():
        model_path = f"{RESULTS}/binary_{embedder}_clustered_model.pt"
        if not os.path.exists(model_path) or (path and not os.path.exists(path)):
            print(f"{embedder}: skipped (model or algae embeddings not there yet)")
            continue

        with open(f"{RESULTS}/binary_{embedder}_clustered_result.json") as f:
            prokaryote = json.load(f)["test_metrics"]
        results[embedder] = {
            **summarise(predict(embedder, ids, model_path), classes, proteomes),
            "prokaryote_specificity": prokaryote["specificity"],
            "prokaryote_accuracy":    prokaryote["accuracy"],
        }

    # ProtT5 also on the unfiltered algae, for the before/after leakage filter comparison
    if "prott5" in results:
        all_ids, all_labels, all_proteome = algae_proteins(filtered=False)
        all_classes = np.array([all_labels[i] for i in all_ids])
        unfiltered = summarise(predict("prott5", all_ids, f"{RESULTS}/binary_prott5_clustered_model.pt"),
                               all_classes, np.array([all_proteome[i] for i in all_ids]))
        results["prott5"]["unfiltered"] = {"n_algae_proteins": len(all_ids), **unfiltered}

    print(f"\n{'Recall on algae':<16}" + "".join(f"{e:>11}" for e in results))
    for c in ["overall"] + CLASSES:
        row = [results[e]["overall" if c == "overall" else "by_class"] for e in results]
        row = [r if c == "overall" else r[c] for r in row]
        print(f"  {c:<14}" + "".join(f"{r['recall']:>11.3f}" for r in row))
    print(f"\n{'Prokaryote':<16}" + "".join(f"{e:>11}" for e in results))
    print(f"  {'specificity':<14}" + "".join(f"{results[e]['prokaryote_specificity']['value']:>11.3f}" for e in results))

    with open(OUT_JSON, "w") as f:
        json.dump({
            "task":     "binary models of each embedder on leakage-filtered eukaryote algae",
            "models":   "binary_<embedder>_clustered_model.pt (clustered binary split, identical for all embedders)",
            "n_algae_proteins": len(ids),
            "results":  results,
        }, f, indent=2)
    print(f"\nResults saved to {OUT_JSON}")
