# Effect of the alkaliphiles whose NaCl optimum also meets the halophile threshold (>= 0.2 M) but that carry only the
# alkaliphile label: five-class results (ProtT5 and composition) and the Asp + Glu share of the alkaliphiles, with and
# without these organisms. Uses the existing models; nothing is retrained.

import csv
import json
from collections import defaultdict

import h5py
import numpy as np
import torch
from sklearn.metrics import f1_score

from acidic_residue_check import properties, summarise
from eval_eukaryote_fiveclass_filtered import CLASS_NAMES, FNN, MODEL_PATH
from eval_fiveclass_by_domain import onehot_model

TEST_FAA     = "dataset_splits/rep_split1_test.faa"
EMBEDDINGS   = "embeddings/split1_30pct.h5"
MASTER_INDEX = "searchable_resource/extremophile_master_index.tsv"
OUT_JSON     = "results/training_results/haloalkaliphile_effect_result.json"

SEED = 42
N_BOOTSTRAP = 1000
HALOPHILE_MIN_NACL = 0.2


def read_test():
    proteins, current = [], None
    with open(TEST_FAA) as f:
        for line in f:
            if line.startswith(">"):
                proteome, description, cls = line[1:].strip().split("||")
                current = {"id": description.split()[0], "proteome": proteome, "class": cls.strip(), "seq": []}
                proteins.append(current)
            else:
                current["seq"].append(line.strip())
    for p in proteins:
        p["seq"] = "".join(p["seq"])
    return proteins


def with_ci(y_true, y_pred):
    def metrics(t, p):
        per_class = f1_score(t, p, average=None, labels=range(5), zero_division=0)
        return {"macro_f1": f1_score(t, p, average="macro"), "f1_alkaliphile": per_class[CLASS_NAMES.index("alkaliphile")]}
    rng = np.random.default_rng(SEED)
    boot = defaultdict(list)
    for _ in range(N_BOOTSTRAP):
        idx = rng.choice(len(y_true), len(y_true), replace=True)
        for k, v in metrics(y_true[idx], y_pred[idx]).items():
            boot[k].append(v)
    return {k: {"value": float(v), "ci_lower": float(np.percentile(boot[k], 2.5)),
                "ci_upper": float(np.percentile(boot[k], 97.5))} for k, v in metrics(y_true, y_pred).items()}


if __name__ == "__main__":
    flagged = {}
    for row in csv.DictReader(open(MASTER_INDEX), delimiter="\t"):
        if (row["domain_scope"] == "prokaryote" and row["phenotype_focus"] == "alkaliphile" and row["nacl_optimum_M"]
                and float(row["nacl_optimum_M"]) >= HALOPHILE_MIN_NACL):
            flagged[row["proteome_id"]] = (row["taxon_name"], float(row["nacl_optimum_M"]))
    print(f"{len(flagged)} alkaliphiles with an NaCl optimum >= {HALOPHILE_MIN_NACL} M:")
    for name, nacl in sorted(flagged.values()):
        print(f"  {name:<36} {nacl:.3f} M")

    proteins = read_test()
    result = {"organisms": {pid: {"taxon": t, "nacl_optimum_M": n} for pid, (t, n) in flagged.items()}, "fiveclass": {}}

    # ProtT5 skips the few test proteins without an embedding, as in eval_fiveclass_by_domain.py
    with h5py.File(EMBEDDINGS, "r") as h5:
        prott5 = [p for p in proteins if p["id"] in h5]
        X_prott5 = np.array([h5[p["id"]][:] for p in prott5], dtype=np.float32)
    model = FNN()
    model.load_state_dict(torch.load(MODEL_PATH, map_location="cpu"))
    composition, aa_composition = onehot_model()
    X_composition = np.array([aa_composition(p["seq"]) for p in proteins], dtype=np.float32)

    print(f"\n{'model':<12} {'test set':<36} {'macro-F1':>22} {'F1 alkaliphile':>22}")
    for name, net, subset, X in (("ProtT5", model, prott5, X_prott5), ("Composition", composition, proteins, X_composition)):
        net.eval()
        with torch.no_grad():
            y_pred = net(torch.from_numpy(X)).argmax(dim=1).numpy()
        y_true = np.array([CLASS_NAMES.index(p["class"]) for p in subset])
        keep = np.array([p["proteome"] not in flagged for p in subset])
        result["fiveclass"][name] = {}
        for label, mask in (("all", np.ones_like(keep)), ("without these organisms", keep)):
            r = with_ci(y_true[mask], y_pred[mask])
            result["fiveclass"][name][label] = {"n_test": int(mask.sum()), **r}
            print(f"{name:<12} {label + f' ({mask.sum():,})':<36} "
                  + "".join(f"{r[k]['value']:>8.3f} [{r[k]['ci_lower']:.3f}, {r[k]['ci_upper']:.3f}]"
                            for k in ("macro_f1", "f1_alkaliphile")))

    groups = {"all alkaliphiles": defaultdict(list), "without these organisms": defaultdict(list),
              "these organisms only": defaultdict(list)}
    for p in proteins:
        props = properties(p["seq"]) if p["class"] == "alkaliphile" else None
        if props is None:
            continue
        groups["all alkaliphiles"][p["proteome"]].append(props)
        groups["these organisms only" if p["proteome"] in flagged else "without these organisms"][p["proteome"]].append(props)
    rng = np.random.default_rng(SEED)
    result["acidic_residues"] = {}
    print(f"\n{'alkaliphile group':<26} {'org.':>4}  {'Asp+Glu %':>20}  {'median pI':>9}")
    for label, per_organism in groups.items():
        s = summarise(per_organism, rng)
        result["acidic_residues"][label] = s
        print(f"{label:<26} {s['n_organisms']:>4}  {s['acidic']['value'] * 100:>5.1f} "
              f"[{s['acidic']['ci_lower'] * 100:.1f}, {s['acidic']['ci_upper'] * 100:.1f}]  {s['pI']['value']:>9.2f}")

    with open(OUT_JSON, "w") as f:
        json.dump(result, f, indent=2)
    print(f"\nSaved to {OUT_JSON}")
