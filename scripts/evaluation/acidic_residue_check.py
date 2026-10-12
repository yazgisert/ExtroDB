# Salt-in check in the sequences: share of acidic (Asp + Glu) and basic (Lys + Arg) residues and the isoelectric point,
# as organism medians per group (halophile subtype x domain, other classes, non-extremophiles, algae) with an
# organism-level bootstrap; also compares bacterial halophile proteins the 5-class model does and does not recognise.

import csv
import json
from collections import defaultdict

import h5py
import numpy as np
import torch
from Bio.SeqUtils.ProtParam import ProteinAnalysis

from eval_eukaryote_fiveclass_filtered import CLASS_NAMES, FNN, MODEL_PATH

TEST_FAA      = "dataset_splits/rep_split1_test.faa"
EMBEDDINGS    = "embeddings/split1_30pct.h5"
DOMAIN_TSV    = "curation/source_tables/proteome_domain_ncbi.tsv"
METADATA      = "dataset_splits/split_input_metadata.tsv"
BINARY_SPLIT  = "dataset_splits/clustered_binary/binary_split.tsv"
NONEXTREM_FAA = "nonextremophile_data/nonextremophile_filtered.faa"
ALGAE_FAA     = "results/homology_inference/filtered_eukaryote_query.faa"
ALGAE_LABELS  = "embeddings/eukaryote_algae_labels.tsv"
MASTER_INDEX  = "searchable_resource/extremophile_master_index.tsv"
OUT_JSON      = "results/training_results/acidic_residue_result.json"

SEED = 42
N_BOOTSTRAP = 1000
STANDARD = set("ACDEFGHIKLMNPQRSTVWY")


def read_faa(path, id_from_header):
    seqs, current, chunks = {}, None, []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if current:
                    seqs[current[0]] = (current[1], "".join(chunks))
                current, chunks = id_from_header(line[1:]), []
            else:
                chunks.append(line)
    if current:
        seqs[current[0]] = (current[1], "".join(chunks))
    return seqs


def properties(seq):
    seq = "".join(aa for aa in seq.upper() if aa in STANDARD)
    if len(seq) < 30:
        return None
    n = len(seq)
    acidic = (seq.count("D") + seq.count("E")) / n
    basic = (seq.count("K") + seq.count("R")) / n
    return {"acidic": acidic, "basic": basic, "net": acidic - basic, "pI": ProteinAnalysis(seq).isoelectric_point()}


def summarise(per_organism, rng):
    # Median per organism, then the mean over organisms; the CI comes from resampling organisms, not proteins
    organisms = sorted(per_organism)
    medians = {k: np.array([np.median([p[k] for p in per_organism[o]]) for o in organisms])
               for k in ("acidic", "basic", "net", "pI")}
    out = {"n_organisms": len(organisms), "n_proteins": sum(len(v) for v in per_organism.values())}
    for k, values in medians.items():
        if len(values) > 1:
            boot = [values[rng.integers(0, len(values), len(values))].mean() for _ in range(N_BOOTSTRAP)]
            ci = (float(np.percentile(boot, 2.5)), float(np.percentile(boot, 97.5)))
        else:
            ci = (None, None)
        out[k] = {"value": float(values.mean()), "ci_lower": ci[0], "ci_upper": ci[1]}
    return out


if __name__ == "__main__":
    rng = np.random.default_rng(SEED)
    domain_of = {r["proteome_id"]: r["domain"] for r in csv.DictReader(open(DOMAIN_TSV), delimiter="\t")}
    subtype_of = {r["proteome_id"]: r["split2_label"] for r in csv.DictReader(open(METADATA), delimiter="\t")}
    taxon_of = {r["proteome_id"]: r["taxon_name"] for r in csv.DictReader(open(MASTER_INDEX), delimiter="\t")}

    groups = defaultdict(lambda: defaultdict(list))

    # Prokaryote extremophiles (5-class test set), with the 5-class prediction for the halophiles
    test = read_faa(TEST_FAA, lambda h: (h.split("||")[1].split()[0], (h.split("||")[0], h.split("||")[2].strip())))
    halophile_ids = []
    for protein_id, ((proteome, cls), seq) in test.items():
        props = properties(seq)
        if props is None:
            continue
        domain = domain_of[proteome]
        if cls == "halophile":
            groups[f"halophile | {subtype_of[proteome]} | {domain}"][proteome].append(props)
            if domain == "bacteria":
                halophile_ids.append((protein_id, props))
        else:
            groups[f"{cls} | {domain}"][proteome].append(props)

    # Non-extremophile test proteins of the binary split (Swiss-Prot bacteria; organism = NCBI taxon in the header)
    with open(BINARY_SPLIT) as f:
        next(f)
        nonext_test = {line.split("\t")[0] for line in f if line.rstrip("\n").endswith("\tnon_extremophile")
                       and line.split("\t")[1] == "test"}

    def nonext_header(h):
        taxid = next((p[3:] for p in h.split() if p.startswith("OX=")), "unknown")
        return h.split()[0], taxid

    for protein_id, (taxid, seq) in read_faa(NONEXTREM_FAA, nonext_header).items():
        if protein_id in nonext_test:
            props = properties(seq)
            if props is not None:
                groups["non-extremophile | Swiss-Prot bacteria"][taxid].append(props)

    # Leakage-filtered algae, per class (organism = proteome)
    algae_label = {}
    with open(ALGAE_LABELS) as f:
        next(f)
        for line in f:
            protein_id, proteome, cls = line.rstrip("\n").split("\t")
            algae_label[protein_id] = (proteome, cls)
    for protein_id, (_, seq) in read_faa(ALGAE_FAA, lambda h: (h.split()[0], None)).items():
        if protein_id in algae_label:
            proteome, cls = algae_label[protein_id]
            props = properties(seq)
            if props is not None:
                groups[f"algae {cls} | eukaryote"][proteome].append(props)

    result = {"groups": {}, "algae_halophile_species": {}, "bacterial_halophiles_by_prediction": {}}
    print(f"{'group':<45} {'org.':>4} {'proteins':>8}  {'Asp+Glu %':>18}  {'Lys+Arg %':>9}  {'median pI':>9}")
    for name in sorted(groups):
        s = summarise(groups[name], rng)
        result["groups"][name] = s
        ci = f"[{s['acidic']['ci_lower'] * 100:.1f}, {s['acidic']['ci_upper'] * 100:.1f}]" if s["acidic"]["ci_lower"] else ""
        print(f"{name:<45} {s['n_organisms']:>4} {s['n_proteins']:>8}  {s['acidic']['value'] * 100:>5.1f} {ci:<12}  "
              f"{s['basic']['value'] * 100:>9.1f}  {s['pI']['value']:>9.2f}")

    print("\nAlgal halophiles per species:")
    for proteome, props in groups["algae halophile | eukaryote"].items():
        acidic, pi = np.median([p["acidic"] for p in props]), np.median([p["pI"] for p in props])
        result["algae_halophile_species"][taxon_of.get(proteome, proteome)] = {"acidic": float(acidic), "pI": float(pi)}
        print(f"  {taxon_of.get(proteome, proteome):<28} Asp+Glu {acidic * 100:.1f}%  median pI {pi:.2f}")

    # Do the bacterial halophile proteins the model gets right look more salt-in-like?
    with h5py.File(EMBEDDINGS, "r") as h5:
        kept = [(i, p) for i, p in halophile_ids if i in h5]
        X = np.array([h5[i][:] for i, _ in kept], dtype=np.float32)
    model = FNN()
    model.load_state_dict(torch.load(MODEL_PATH, map_location="cpu"))
    model.eval()
    with torch.no_grad():
        called_halophile = model(torch.from_numpy(X)).argmax(dim=1).numpy() == CLASS_NAMES.index("halophile")
    print("\nBacterial halophile test proteins by 5-class prediction:")
    for label, mask in (("predicted halophile", called_halophile), ("predicted other class", ~called_halophile)):
        props = [p for (_, p), m in zip(kept, mask) if m]
        entry = {k: float(np.median([p[k] for p in props])) for k in ("acidic", "basic", "net", "pI")}
        entry["n_proteins"] = len(props)
        result["bacterial_halophiles_by_prediction"][label] = entry
        print(f"  {label:<22} n={len(props):>5}  Asp+Glu {entry['acidic'] * 100:.1f}%  Lys+Arg {entry['basic'] * 100:.1f}%  "
              f"median pI {entry['pI']:.2f}")

    with open(OUT_JSON, "w") as f:
        json.dump(result, f, indent=2)
    print(f"\nResults saved to {OUT_JSON}")
