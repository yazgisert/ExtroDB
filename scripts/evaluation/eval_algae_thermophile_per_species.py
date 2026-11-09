# Thermophily per algal species: THERMO's predicted probability (one prediction per proteome, from
# run_thermo_eukaryote_eval.py) next to the share of the species' proteins that the ProtT5 five-class and six-class
# models assign to the thermophile class. Uses the leakage-filtered algae proteins with a ProtT5 embedding, as all other
# algae results. The random baseline is the share of thermophile labels in each model's training set, with a 95%
# bootstrap interval for the number of proteins of the species.
# Run from the repository root (biocentral env): python scripts/evaluation/eval_algae_thermophile_per_species.py

import json
import sys
from collections import Counter
from pathlib import Path

import numpy as np
import torch

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "training"))
import eval_eukaryote_embedders as algae_eval
import eval_eukaryote_fiveclass_filtered as five
import train_sixclass_extremophile as six

RESULTS = Path("results/training_results")
OUT_JSON = RESULTS / "algae_thermophile_per_species_result.json"
SEED, N_BOOTSTRAP = 42, 1000


def thermophile_predictions(model, X, class_names):
    model.eval()
    with torch.no_grad():
        return model(torch.from_numpy(X)).argmax(dim=1).numpy() == class_names.index("thermophile")


def random_share_ci(prior, n, rng):
    # Share of n proteins labelled thermophile when each label is drawn with the training probability
    draws = rng.binomial(n, prior, N_BOOTSTRAP) / n
    return float(np.percentile(draws, 2.5)), float(np.percentile(draws, 97.5))


if __name__ == "__main__":
    ids, labels, proteome = algae_eval.algae_proteins()
    X = algae_eval.features("prott5", ids)

    model5 = five.FNN()
    model5.load_state_dict(torch.load(five.MODEL_PATH, map_location="cpu"))
    model6 = six.FNN()
    model6.load_state_dict(torch.load(six.MODEL_PATH, map_location="cpu"))
    pred5 = thermophile_predictions(model5, X, five.CLASS_NAMES)
    pred6 = thermophile_predictions(model6, X, six.CLASS_NAMES)

    # Training share of the thermophile class: five-class Split 1 (balanced) and six-class clustered split
    prior5 = 1 / len(five.CLASS_NAMES)
    six_train = Counter(line.rstrip("\n").split("\t")[2] for line in list(open(six.SPLIT_TSV))[1:]
                        if line.split("\t")[1] == "train")
    prior6 = six_train["thermophile"] / sum(six_train.values())

    thermo = {p["proteome_id"]: p for p in json.load(open(RESULTS / "thermo_eukaryote_eval_result.json"))["per_proteome"]}
    rng = np.random.default_rng(SEED)
    species = []
    proteomes = np.array([proteome[i] for i in ids])
    for pid, t in thermo.items():
        mask = proteomes == pid
        n = int(mask.sum())
        species.append({"proteome_id": pid, "species": t["species"], "broad_class": t["broad_class"],
                        "thermophile": bool(t["true_label"]), "n_proteins": n,
                        "thermo_probability": t["proba_thermophile"],
                        "fiveclass_thermophile_share": float(pred5[mask].mean()),
                        "sixclass_thermophile_share": float(pred6[mask].mean()),
                        "random_fiveclass_ci": random_share_ci(prior5, n, rng),
                        "random_sixclass_ci": random_share_ci(prior6, n, rng)})
        s = species[-1]
        print(f"{s['species']:<32} {s['broad_class']:<13} n={n:>6,}  THERMO {s['thermo_probability']:.3f}  "
              f"5-class {s['fiveclass_thermophile_share']:.3f}  6-class {s['sixclass_thermophile_share']:.3f}")

    result = {"random_fiveclass": prior5, "random_sixclass": prior6, "species": species}
    for key in ("thermo_probability", "fiveclass_thermophile_share", "sixclass_thermophile_share"):
        hot = [s[key] for s in species if s["thermophile"]]
        rest = [s[key] for s in species if not s["thermophile"]]
        ranks = sorted(species, key=lambda s: -s[key])
        result[key] = {"mean_thermophilic": float(np.mean(hot)), "mean_other": float(np.mean(rest)),
                       "thermophilic_in_top4": sum(s["thermophile"] for s in ranks[:4])}
        print(f"{key:<30} thermophilic algae {np.mean(hot):.3f}, other algae {np.mean(rest):.3f}, "
              f"thermophilic among top 4: {result[key]['thermophilic_in_top4']}")
    print(f"random: five-class {prior5:.3f}, six-class {prior6:.3f}")

    with open(OUT_JSON, "w") as f:
        json.dump(result, f, indent=2)
    print(f"\nSaved to {OUT_JSON}")
