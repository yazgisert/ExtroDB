# Checks how many non-extremophile proteins of the clustered binary split come from thermophilic (OGT >= 45 C) or
# psychrophilic (OGT <= 15 C) species according to TEMPURA, and how the binary ProtT5 model performs on the test set
# without them. Species are matched by genus and species name; for species with several strains in TEMPURA, the
# optimal growth temperatures are averaged.
# Run from the repository root: python scripts/evaluation/check_nonextremophile_contamination.py

import csv
import json
import re
from collections import Counter, defaultdict

import h5py
import numpy as np
import torch
import torch.nn as nn
from sklearn.metrics import matthews_corrcoef

SPLIT_TSV     = "dataset_splits/clustered_binary/binary_split.tsv"
NONEXTREM_FAA = "nonextremophile_data/nonextremophile_filtered.faa"
TEMPURA_CSV   = "curation/source_tables/200617_TEMPURA.csv"
EXTREM_H5     = "embeddings/split1_30pct.h5"
NONEXTREM_H5  = "embeddings/nonextremophile.h5"
MODEL_PATH    = "results/training_results/binary_prott5_clustered_model.pt"
OUT_JSON      = "results/training_results/nonextremophile_contamination_result.json"

THERMOPHILE_MIN  = 45.0
PSYCHROPHILE_MAX = 15.0


class FNN(nn.Module):
    def __init__(self, n_features=1024, n_classes=2, dropout_rate=0.25):
        super().__init__()
        self.fc1     = nn.Linear(n_features, 32)
        self.dropout = nn.Dropout(dropout_rate)
        self.fc2     = nn.Linear(32, n_classes)

    def forward(self, x):
        return self.fc2(self.dropout(torch.relu(self.fc1(x))))


def species_name(organism):
    # "Thermus thermophilus (strain HB8)" -> "Thermus thermophilus"
    return " ".join(re.sub(r"\(.*?\)", "", organism).replace('"', "").split()[:2])


def read_split():
    extrem, nonextrem = defaultdict(list), defaultdict(list)
    with open(SPLIT_TSV) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            (nonextrem if row["label"] == "non_extremophile" else extrem)[row["split"]].append(row["protein_id"])
    return extrem, nonextrem


def read_species():
    species = {}
    with open(NONEXTREM_FAA) as f:
        for line in f:
            if line.startswith(">"):
                match = re.search(r" OS=(.+?) OX=", line)
                species[line[1:].split()[0]] = species_name(match.group(1)) if match else ""
    return species


def read_tempura():
    temperatures = defaultdict(list)
    with open(TEMPURA_CSV) as f:
        for row in csv.DictReader(f):
            try:
                temperatures[species_name(row["genus_and_species"])].append(float(row["Topt_ave"]))
            except ValueError:
                continue
    return {name: float(np.mean(values)) for name, values in temperatures.items()}


def embeddings(ids, path):
    with h5py.File(path, "r") as f:
        return np.array([f[i][:] for i in ids], dtype=np.float32)


if __name__ == "__main__":
    extrem, nonextrem = read_split()
    species = read_species()
    tempura = read_tempura()

    def phenotype(protein_id):
        ogt = tempura.get(species[protein_id])
        if ogt is None:
            return None
        if ogt >= THERMOPHILE_MIN:
            return "thermophile"
        if ogt <= PSYCHROPHILE_MAX:
            return "psychrophile"
        return "other"

    print("Non-extremophile proteins from thermophilic or psychrophilic species (TEMPURA)")
    per_split, all_flagged, n_all, n_in_tempura = {}, [], 0, 0
    for split in ("train", "val", "test"):
        ids = nonextrem[split]
        labels = [phenotype(p) for p in ids]
        flagged = [p for p, label in zip(ids, labels) if label in ("thermophile", "psychrophile")]
        per_split[split] = {"n": len(ids), "flagged": len(flagged),
                            "thermophile": labels.count("thermophile"), "psychrophile": labels.count("psychrophile")}
        all_flagged += flagged
        n_all += len(ids)
        n_in_tempura += sum(label is not None for label in labels)
        print(f"  {split:<5} {len(flagged):>4} / {len(ids):>5}  ({100 * len(flagged) / len(ids):.1f}%)")
    print(f"  total {len(all_flagged):>4} / {n_all:>5}  ({100 * len(all_flagged) / n_all:.1f}%); "
          f"{n_in_tempura:,} proteins have a species listed in TEMPURA")

    top_species = Counter(species[p] for p in all_flagged).most_common(10)
    print("\nMost frequent species:")
    for name, count in top_species:
        print(f"  {name:<40} {count:>4}  (OGT {tempura[name]:.0f} C)")

    # Binary ProtT5 model on the test set, with and without the flagged proteins
    model = FNN()
    model.load_state_dict(torch.load(MODEL_PATH, map_location="cpu"))
    model.eval()
    ext_ids, non_ids = extrem["test"], nonextrem["test"]
    X = np.vstack([embeddings(ext_ids, EXTREM_H5), embeddings(non_ids, NONEXTREM_H5)])
    y = np.array([1] * len(ext_ids) + [0] * len(non_ids))
    with torch.no_grad():
        y_pred = model(torch.from_numpy(X)).argmax(dim=1).numpy()

    flagged_test = set(all_flagged) & set(non_ids)
    is_flagged = np.array([False] * len(ext_ids) + [p in flagged_test for p in non_ids])
    is_non = y == 0
    keep = ~is_flagged
    effect = {
        "mcc_all":                         float(matthews_corrcoef(y, y_pred)),
        "mcc_without_flagged":             float(matthews_corrcoef(y[keep], y_pred[keep])),
        "specificity_all":                 float(1 - y_pred[is_non].mean()),
        "specificity_without_flagged":     float(1 - y_pred[is_non & keep].mean()),
        "predicted_extremophile_flagged":  float(y_pred[is_flagged].mean()),
        "predicted_extremophile_other_non": float(y_pred[is_non & keep].mean()),
    }
    print(f"\nBinary ProtT5 model on the test set ({len(flagged_test)} flagged proteins)")
    print(f"  MCC          all {effect['mcc_all']:.3f}   without flagged {effect['mcc_without_flagged']:.3f}")
    print(f"  specificity  all {effect['specificity_all']:.3f}   without flagged {effect['specificity_without_flagged']:.3f}")
    print(f"  predicted as extremophile: flagged {effect['predicted_extremophile_flagged']:.3f}   "
          f"other non-extremophiles {effect['predicted_extremophile_other_non']:.3f}")

    result = {
        "task":       "non-extremophile proteins of the clustered binary split from thermophiles or psychrophiles",
        "source":     TEMPURA_CSV,
        "thresholds": {"thermophile_min_ogt": THERMOPHILE_MIN, "psychrophile_max_ogt": PSYCHROPHILE_MAX},
        "n_nonextremophile": n_all,
        "n_with_tempura_species": n_in_tempura,
        "n_flagged":  len(all_flagged),
        "per_split":  per_split,
        "top_species": [{"species": name, "n_proteins": count, "ogt": tempura[name]} for name, count in top_species],
        "binary_prott5_test": effect,
    }
    with open(OUT_JSON, "w") as f:
        json.dump(result, f, indent=2)
    print(f"\nSaved to {OUT_JSON}")
