# Archaea bias check for the clustered binary models: the non-extremophiles are bacteria only (Swiss-Prot
# taxonomy 2), while a quarter of the extremophile organisms are archaea. Evaluates every binary model on the full
# test set and on a bacteria-only test set (archaeal extremophiles removed), with 1000 bootstrap CIs, and reports
# the ProtT5 recall per class and domain. Domains come from NCBI Taxonomy (fetched once and cached).

import csv
import json
import sys
import urllib.request
import xml.etree.ElementTree as ET
from collections import defaultdict
from pathlib import Path

import h5py
import numpy as np
import torch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "training"))
from train_binary_embedder import (EMBEDDINGS, EXTREM_FAA, FNN, NONEXTREM_FAA, amino_acid_composition, read_faa,
                                   read_split, split_id, with_ci)

MASTER_INDEX = "searchable_resource/extremophile_master_index.tsv"
DOMAIN_TSV   = "curation/source_tables/proteome_domain_ncbi.tsv"
OUT_JSON     = "results/training_results/binary_domain_bias_result.json"


def proteome_domains():
    if not Path(DOMAIN_TSV).exists():
        rows = [r for r in csv.DictReader(open(MASTER_INDEX), delimiter="\t") if r["domain_scope"] == "prokaryote"]
        taxids = sorted({r["ncbi_taxid"].strip() for r in rows if r["ncbi_taxid"].strip()})
        url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=taxonomy&id=" + ",".join(taxids)
        domain_of_taxid = {}
        for taxon in ET.fromstring(urllib.request.urlopen(url, timeout=120).read()).findall("Taxon"):
            lineage = taxon.findtext("Lineage") or ""
            domain = "archaea" if "Archaea" in lineage else "bacteria" if "Bacteria" in lineage else "unknown"
            for taxid in [taxon.findtext("TaxId")] + [t.text for t in taxon.findall("AkaTaxIds/TaxId")]:
                domain_of_taxid[taxid] = domain
        with open(DOMAIN_TSV, "w") as f:
            f.write("proteome_id\ttaxon_name\tncbi_taxid\tdomain\n")
            for r in rows:
                taxid = r["ncbi_taxid"].strip()
                f.write(f"{r['proteome_id']}\t{r['taxon_name']}\t{taxid}\t{domain_of_taxid.get(taxid, 'unknown')}\n")
        print(f"  NCBI domains saved to {DOMAIN_TSV}")
    return {r["proteome_id"]: r["domain"] for r in csv.DictReader(open(DOMAIN_TSV), delimiter="\t")}


def protein_info(domains):
    info = {}
    for path in EXTREM_FAA:
        with open(path) as f:
            for line in f:
                if line.startswith(">"):
                    parts = line[1:].strip().split("||")
                    info[parts[1].split()[0]] = (domains[parts[0]], parts[2])
    return info


def features(embedder, extrem_ids, nonextrem_ids):
    if embedder == "onehot":
        seqs = {}
        for path in EXTREM_FAA:
            seqs.update(read_faa(path, split_id))
        seqs.update(read_faa(NONEXTREM_FAA, lambda h: h.split()[0]))
        return np.array([amino_acid_composition(seqs[i]) for i in extrem_ids + nonextrem_ids], dtype=np.float32)
    extrem_h5, nonextrem_h5 = EMBEDDINGS[embedder]
    with h5py.File(extrem_h5, "r") as f:
        X = [f[i][:] for i in extrem_ids]
    with h5py.File(nonextrem_h5, "r") as f:
        X += [f[i][:] for i in nonextrem_ids]
    return np.array(X, dtype=np.float32)


if __name__ == "__main__":
    domains = proteome_domains()
    info = protein_info(domains)
    extrem, nonextrem = read_split()
    extrem_ids, nonextrem_ids = extrem["test"], nonextrem["test"]
    is_archaeal = np.array([info[i][0] == "archaea" for i in extrem_ids] + [False] * len(nonextrem_ids))
    y_true = np.array([1] * len(extrem_ids) + [0] * len(nonextrem_ids))
    print(f"Test set: {len(extrem_ids):,} extremophiles ({is_archaeal.sum():,} archaeal), "
          f"{len(nonextrem_ids):,} non-extremophiles (all bacterial)")

    result = {"n_test_extremophile": len(extrem_ids), "n_test_archaeal_extremophile": int(is_archaeal.sum()),
              "n_test_nonextremophile": len(nonextrem_ids), "models": {}}
    for embedder in EMBEDDINGS:
        X = features(embedder, extrem_ids, nonextrem_ids)
        model = FNN(n_features=X.shape[1])
        model.load_state_dict(torch.load(f"results/training_results/binary_{embedder}_clustered_model.pt",
                                         map_location="cpu"))
        model.eval()
        with torch.no_grad():
            y_pred = model(torch.from_numpy(X)).argmax(dim=1).numpy()

        bacterial = ~is_archaeal
        archaeal_extrem = is_archaeal & (y_true == 1)
        result["models"][embedder] = {
            "all":                    with_ci(y_true, y_pred),
            "bacteria_only":          with_ci(y_true[bacterial], y_pred[bacterial]),
            "recall_archaeal_extrem": float(y_pred[archaeal_extrem].mean()),
        }
        if embedder == "prott5":
            per_class = defaultdict(lambda: defaultdict(list))
            for i, pid in enumerate(extrem_ids):
                domain, cls = info[pid]
                per_class[cls][domain].append(y_pred[i])
            result["prott5_recall_per_class"] = {
                cls: {d: {"n": len(p), "recall": float(np.mean(p))} for d, p in by_domain.items()}
                for cls, by_domain in per_class.items()
            }

        m = result["models"][embedder]
        print(f"\n{embedder}")
        for subset in ("all", "bacteria_only"):
            s = m[subset]
            print(f"  {subset:<14} MCC {s['mcc']['value']:.3f} [{s['mcc']['ci_lower']:.3f}, {s['mcc']['ci_upper']:.3f}]"
                  f"  accuracy {s['accuracy']['value']:.3f}  recall {s['recall']['value']:.3f}"
                  f"  specificity {s['specificity']['value']:.3f}")
        print(f"  recall on archaeal extremophiles: {m['recall_archaeal_extrem']:.3f}")

    with open(OUT_JSON, "w") as f:
        json.dump(result, f, indent=2)
    print(f"\nResults saved to {OUT_JSON}")
