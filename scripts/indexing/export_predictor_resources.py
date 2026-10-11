# Exports the predictor inputs for the website: one gzipped FASTA per model and split (5-/12-class, clustered binary
# and 6-class, filtered algae) plus a README for Zenodo, and manifest.json / results_summary.json for predictors.html.

import csv
import gzip
import json
from collections import Counter, defaultdict
from pathlib import Path

import h5py

OUT_DIR    = Path("predictor_resources")
WEB_DIR    = Path("searchable_resource_app/predictors")
ZENODO_DIR = OUT_DIR / "zenodo"
RESULTS    = Path("results/training_results")

FIVECLASS_SPLITS = {
    "fiveclass_split1_30pct": ("rep_split1",       "5 broad classes, 30% identity clustering"),
    "fiveclass_split1_20pct": ("rep_split1_20pct", "5 broad classes, 20% identity clustering"),
    "twelveclass_split2_30pct": ("rep_split2",       "12 subtypes, 30% identity clustering"),
    "twelveclass_split2_20pct": ("rep_split2_20pct", "12 subtypes, 20% identity clustering"),
}
RESULT_PREFIX = {
    "fiveclass_split1_30pct": "split1_30pct", "fiveclass_split1_20pct": "split1_20pct",
    "twelveclass_split2_30pct": "split2_30pct", "twelveclass_split2_20pct": "split2_20pct",
}
NONEXTREM_FAA = "nonextremophile_data/nonextremophile_filtered.faa"
ALGAE_FAA     = "results/homology_inference/filtered_eukaryote_query.faa"
ALGAE_H5      = "embeddings/eukaryote_algae.h5"
ALGAE_LABELS  = "embeddings/eukaryote_algae_labels.tsv"
BLACKLIST     = "results/eukaryote_leakage/flagged_protein_ids.txt"
MASTER_INDEX  = "searchable_resource/extremophile_master_index.tsv"
EMBEDDERS     = {"prott5": "ProtT5", "esm2_650m": "ESM2 650M", "esm2_8m": "ESM2 8M", "onehot": "One-hot"}
ALGAE_CLASSES = ["acidophile", "halophile", "psychrophile", "thermophile"]


def read_faa(path, parse_header):
    records = {}
    current, seq = None, []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if current:
                    records[current[0]] = (current[1], "".join(seq))
                current, seq = parse_header(line[1:]), []
            else:
                seq.append(line)
    if current:
        records[current[0]] = (current[1], "".join(seq))
    return records


def split_header(header):
    proteome, rest = header.split("||", 1)
    return rest.split()[0], {"label": header.rsplit("||", 1)[1].strip(), "proteome": proteome}


def write_fasta(name, entries):
    path = ZENODO_DIR / f"{name}.fasta.gz"
    with gzip.open(path, "wt") as f:
        for protein_id, info, seq in entries:
            fields = " ".join(f"{k}={v}" for k, v in info.items())
            f.write(f">{protein_id} {fields}\n{seq}\n")
    counts = Counter(info["label"] for _, info, _ in entries)
    return {"file": path.name, "n_proteins": len(entries), "labels": dict(sorted(counts.items())),
            "size_mb": round(path.stat().st_size / 1e6, 2)}


def fiveclass_files():
    files = {}
    for name, (prefix, description) in FIVECLASS_SPLITS.items():
        files[name] = {"description": description, "splits": {}}
        for split in ("train", "val", "test"):
            records = read_faa(f"dataset_splits/{prefix}_{split}.faa", split_header)
            entries = [(pid, info, seq) for pid, (info, seq) in sorted(records.items())]
            files[name]["splits"][split] = write_fasta(f"{name}_{split}", entries)
    return files


def clustered_files(extremophile_records):
    nonext = read_faa(NONEXTREM_FAA, lambda h: (h.split()[0], {}))
    files = {}
    for name, tsv, description in [
        ("binary_clustered", "dataset_splits/clustered_binary/binary_split.tsv",
         "extremophile vs. non-extremophile; 13,031 non-extremophile cluster representatives (30%) + as many extremophiles, equal per class"),
        ("sixclass_clustered", "dataset_splits/clustered_binary/sixclass_split.tsv",
         "5 extremophile classes + non-extremophile, all classes equal size; test set identical to binary_clustered"),
    ]:
        rows = defaultdict(list)
        with open(tsv) as f:
            next(f)
            for line in f:
                protein_id, split, label = line.rstrip("\n").split("\t")
                rows[split].append((protein_id, label))
        files[name] = {"description": description, "splits": {}}
        for split in ("train", "val", "test"):
            entries = []
            for protein_id, label in sorted(rows[split]):
                if label == "non_extremophile":
                    seq = nonext[protein_id][1]
                    info = {"label": label, "binary": "non_extremophile", "source": "UniProt_SwissProt"}
                else:
                    meta, seq = extremophile_records[protein_id]
                    info = {"label": label, "binary": "extremophile", "proteome": meta["proteome"]}
                entries.append((protein_id, info, seq))
            files[name]["splits"][split] = write_fasta(f"{name}_{split}", entries)
    return files


def algae_file():
    with open(BLACKLIST) as f:
        blacklist = {line.strip() for line in f if line.strip()}
    with h5py.File(ALGAE_H5, "r") as f:
        embedded = set(f.keys())
    with open(MASTER_INDEX) as f:
        species = {row["proteome_id"]: row["taxon_name"] for row in csv.DictReader(f, delimiter="\t")}
    seqs = read_faa(ALGAE_FAA, lambda h: (h.split()[0], {}))

    entries = []
    with open(ALGAE_LABELS) as f:
        next(f)
        for line in f:
            protein_id, proteome, label = line.rstrip("\n").split("\t")
            if protein_id in embedded and protein_id not in blacklist and label in ALGAE_CLASSES:
                info = {"label": label, "proteome": proteome, "species": species.get(proteome, proteome).replace(" ", "_")}
                entries.append((protein_id, info, seqs[protein_id][1]))
    entries.sort()
    return {"algae_filtered": {
        "description": "eukaryotic algae, test only; proteins with >= 30% identity to prokaryote train + val removed",
        "splits": {"test": write_fasta("algae_filtered_test", entries)},
    }}


def metric(m, key):
    if "value" in m[key]:
        return {"value": m[key]["value"], "ci_lower": m[key]["ci_lower"], "ci_upper": m[key]["ci_upper"]}
    return {"value": m[key], "ci_lower": None, "ci_upper": None}


def results_summary():
    summary = {"fiveclass": {}, "binary_clustered": {}, "binary_bacteria_only": {}, "sixclass_clustered": {}, "algae": {}, "homology": {}}

    suffix = {"prott5": "", "esm2_650m": "_esm2", "esm2_8m": "_esm2_8m", "onehot": "_onehot"}
    for name, prefix in RESULT_PREFIX.items():
        summary["fiveclass"][name] = {}
        for key, label in EMBEDDERS.items():
            r = json.load(open(RESULTS / f"{prefix}{suffix[key]}_result.json"))
            m, b = r["test_metrics"], r["test_metrics_bootstrap"]
            summary["fiveclass"][name][label] = {
                k: {"value": m[src], "ci_lower": b[src]["lower"], "ci_upper": b[src]["upper"]}
                for k, src in (("accuracy", "accuracy"), ("macro_f1", "macro-f1_score"), ("mcc", "matthews-corr-coeff"))
            } | {"n_test": r["n_test"]}

    for key, label in EMBEDDERS.items():
        r = json.load(open(RESULTS / f"binary_{key}_clustered_result.json"))
        summary["binary_clustered"][label] = {k: metric(r["test_metrics"], k)
                                              for k in ("accuracy", "specificity", "recall", "f1_extremophile", "mcc")}

    # Non-extremophiles are bacteria only, so the same models are also scored without the archaeal extremophiles
    bias = json.load(open(RESULTS / "binary_domain_bias_result.json"))["models"]
    summary["binary_bacteria_only"] = {label: {k: metric(bias[key]["bacteria_only"], k)
                                               for k in ("accuracy", "specificity", "recall", "mcc")}
                                       for key, label in EMBEDDERS.items()}

    six = json.load(open(RESULTS / "sixclass_clustered_result.json"))
    five = json.load(open(RESULTS / "fiveclass_binarized_clustered_result.json"))
    summary["sixclass_clustered"] = {
        "Binarized 5-class": {k: metric(five["binarized_fiveclass"], k) for k in ("accuracy", "specificity", "recall", "mcc")},
        "Binarized 6-class": {k: metric(six["binarized_sixclass"], k) for k in ("accuracy", "specificity", "recall", "mcc")},
        "Dedicated binary":  {k: metric(six["dedicated_binary"], k) for k in ("accuracy", "specificity", "recall", "mcc")},
    }

    algae = json.load(open(RESULTS / "eukaryote_embedder_clustered_result.json"))["results"]
    for key, label in EMBEDDERS.items():
        r = algae[key]
        summary["algae"][label] = {"overall": r["overall"]["recall"],
                                   **{c: r["by_class"][c]["recall"] for c in ALGAE_CLASSES}}

    hom5 = json.load(open(RESULTS / "homology_fiveclass_result.json"))
    hom2 = json.load(open(RESULTS / "homology_binary_clustered_result.json"))
    summary["homology"] = {
        "fiveclass_split1_30pct": {"homology_accuracy": hom5["homology"]["accuracy"]["value"],
                                   "prott5_accuracy": hom5["prott5"]["accuracy"]["value"],
                                   "prott5_accuracy_without_hit": hom5["subsets"]["without_hit"]["prott5_accuracy"]},
        "binary_clustered": {"homology_accuracy": hom2["homology"]["accuracy"]["value"],
                             "prott5_accuracy": hom2["prott5"]["accuracy"]["value"],
                             "prott5_accuracy_without_hit": hom2["subsets"]["without_hit"]["prott5_accuracy"]},
    }
    return summary


README = """# Predictor input FASTAs

Input data for the extremophile predictors of the thesis. Every model was trained
and evaluated with ProtT5, ESM2 650M, ESM2 8M and one-hot (amino acid composition)
embeddings of exactly these sequences.

Header format: `>protein_id label=<class> key=value ...`

| Files | Model | Content |
|---|---|---|
{rows}

Notes
- Extremophile proteins are MMseqs2 cluster representatives (one per cluster), split
  70/10/20 per class and downsampled to the smallest class.
- 60 proteins in the 5-class files have no ProtT5 embedding; the ProtT5 models were
  trained without them.
- `algae_filtered_test` contains only proteins without >= 30% identity to any
  prokaryote training or validation protein.
"""


if __name__ == "__main__":
    ZENODO_DIR.mkdir(parents=True, exist_ok=True)

    extremophile_records = {}
    for split in ("train", "val", "test"):
        extremophile_records.update(read_faa(f"dataset_splits/rep_split1_{split}.faa", split_header))

    print("Writing FASTAs ...")
    manifest = {**fiveclass_files(), **clustered_files(extremophile_records), **algae_file()}

    rows = []
    for name, model in manifest.items():
        for split, info in model["splits"].items():
            rows.append(f"| `{info['file']}` | {name} ({split}) | {info['n_proteins']:,} proteins · {model['description']} |")
            print(f"  {info['file']:<42} {info['n_proteins']:>9,}  {info['size_mb']:>6} MB")
    (ZENODO_DIR / "README.md").write_text(README.format(rows="\n".join(rows)))

    # Same two files go next to predictors.html so the page can load them
    WEB_DIR.mkdir(parents=True, exist_ok=True)
    summary = results_summary()
    for folder in (OUT_DIR, WEB_DIR):
        with open(folder / "manifest.json", "w") as f:
            json.dump(manifest, f, indent=2)
        with open(folder / "results_summary.json", "w") as f:
            json.dump(summary, f, indent=2)

    total = sum(p.stat().st_size for p in ZENODO_DIR.glob("*.fasta.gz")) / 1e6
    print(f"\n{len(rows)} FASTA files, {total:.0f} MB in total to {ZENODO_DIR}")
    print(f"Manifest and results summary to {OUT_DIR} and {WEB_DIR}")
