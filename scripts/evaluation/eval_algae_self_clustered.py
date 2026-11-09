# Clusters the leakage-filtered algae proteins against each other at 30% identity (MMseqs2 easy-cluster, as for the
# prokaryote splits) and evaluates the binary models (all embedders) and the 5-class ProtT5 model on the cluster
# representatives only, to check whether redundancy within the algae inflates the algae results.

import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import torch

sys.path.insert(0, str(Path(__file__).parent))
import eval_eukaryote_embedders as emb

FIVECLASS_MODEL = "results/training_results/fiveclass_extremophile_model.pt"
FIVECLASS_NAMES = ["acidophile", "alkaliphile", "halophile", "psychrophile", "thermophile"]
OUT_DIR  = Path("results/eukaryote_leakage/algae_self_cluster_30pct")
OUT_JSON = "results/training_results/algae_self_clustered_result.json"
MIN_SEQ_ID = 0.3


def cluster(ids):
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    rep_faa = OUT_DIR / "algae_rep_seq.fasta"
    if not rep_faa.exists():
        seqs = emb.read_faa(emb.ALGAE_FAA)
        faa = OUT_DIR / "algae_filtered.faa"
        with open(faa, "w") as f:
            for i in ids:
                f.write(f">{i}\n{seqs[i]}\n")
        subprocess.run(["mmseqs", "easy-cluster", str(faa), str(OUT_DIR / "algae"), str(OUT_DIR / "tmp"),
                        "--min-seq-id", str(MIN_SEQ_ID), "--threads", "8"], check=True, stdout=subprocess.DEVNULL)
    return {line[1:].split()[0] for line in open(rep_faa) if line.startswith(">")}


def fiveclass_correct(ids, labels):
    X = emb.features("prott5", ids)
    model = emb.FNN(n_features=X.shape[1], n_classes=5)
    model.load_state_dict(torch.load(FIVECLASS_MODEL, map_location="cpu"))
    model.eval()
    with torch.no_grad():
        pred = model(torch.from_numpy(X)).argmax(dim=1).numpy()
    correct = (pred == np.array([FIVECLASS_NAMES.index(labels[i]) for i in ids])).astype(int)
    return emb.recall_with_ci(correct)


if __name__ == "__main__":
    ids, labels, proteome = emb.algae_proteins()
    reps = cluster(ids)
    rep_ids = [i for i in ids if i in reps]
    print(f"{len(ids):,} filtered algae proteins, {len(rep_ids):,} cluster representatives at 30%\n")

    result = {"n_all": len(ids), "n_representatives": len(rep_ids), "min_seq_id": MIN_SEQ_ID, "binary": {}}
    for embedder in emb.ALGAE_EMBEDDINGS:
        model_path = f"{emb.RESULTS}/binary_{embedder}_clustered_model.pt"
        entry = {}
        for name, subset in [("all", ids), ("representatives", rep_ids)]:
            pred = (emb.predict(embedder, subset, model_path) >= 0.5).astype(int)
            classes = np.array([labels[i] for i in subset])
            entry[name] = {"overall": emb.recall_with_ci(pred),
                           "by_class": {c: emb.recall_with_ci(pred[classes == c]) for c in emb.CLASSES}}
        result["binary"][embedder] = entry
        print(f"binary {embedder:<10} recall all {entry['all']['overall']['recall']:.3f}   "
              f"representatives {entry['representatives']['overall']['recall']:.3f}")
        for c in emb.CLASSES:
            print(f"    {c:<13} {entry['all']['by_class'][c]['recall']:.3f}   "
                  f"{entry['representatives']['by_class'][c]['recall']:.3f}")

    result["fiveclass_prott5"] = {"all": fiveclass_correct(ids, labels),
                                  "representatives": fiveclass_correct(rep_ids, labels)}
    print(f"\n5-class ProtT5 correct class: all {result['fiveclass_prott5']['all']['recall']:.3f}   "
          f"representatives {result['fiveclass_prott5']['representatives']['recall']:.3f}")

    with open(OUT_JSON, "w") as f:
        json.dump(result, f, indent=2)
    print(f"\nSaved to {OUT_JSON}")
