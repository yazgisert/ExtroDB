# Prepares biotrainer runs for every predictor: one input FASTA in biotrainer format (>ID TARGET=<class> SET=<split>),
# the existing embeddings converted to biotrainer's h5 format (one dataset per unique sequence, keyed by the sequence
# hash biotrainer uses, with an "original_id" attribute), and one config.yml per embedder with the same FNN setup as
# before. Proteins missing from any embedding file are left out, so all embedders are trained on the same proteins.
# Run from the repository root (biocentral env): python scripts/biotrainer/prepare_biotrainer_runs.py

import hashlib
import sys
from pathlib import Path

import h5py
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "training"))
from train_onehot import aa_composition

RUNS_DIR      = Path("biotrainer_runs")
EMBED_DIR     = RUNS_DIR / "embeddings"
SPLITS        = Path("dataset_splits")
NONEXTREM_FAA = Path("nonextremophile_data/nonextremophile_filtered.faa")

EMBEDDERS = {
    "prott5":      {"suffix": "",         "name": "Rostlab/prot_t5_xl_uniref50"},
    "esm2_650m":   {"suffix": "_esm2",    "name": "facebook/esm2_t33_650M_UR50D"},
    "esm2_8m":     {"suffix": "_esm2_8m", "name": "facebook/esm2_t6_8M_UR50D"},
    "composition": {"suffix": None,       "name": "amino acid composition (20 relative frequencies)"},
}

SPLIT_PREDICTORS = {
    "fiveclass_split1_30pct":   ("rep_split1",       "split1_30pct"),
    "fiveclass_split1_20pct":   ("rep_split1_20pct", "split1_20pct"),
    "twelveclass_split2_30pct": ("rep_split2",       "split2_30pct"),
    "twelveclass_split2_20pct": ("rep_split2_20pct", "split2_20pct"),
}

CONFIG = """input_file: ../input.fasta
embeddings_file: ../../embeddings/{predictor}_{embedder}.h5
protocol: sequence_to_class
model_choice: FNN
optimizer_choice: adam
loss_choice: cross_entropy_loss
learning_rate: 1e-3
batch_size: 128
num_epochs: 200
patience: 10
epsilon: 1e-3
dropout_rate: 0.25
use_class_weights: False
seed: 42
device: cpu
bootstrapping_iterations: 1000
external_writer: none
disable_pytorch_compile: True
# Embeddings: {embedder_name}
"""


def read_faa(path, parse_id):
    seqs, current, chunks = {}, None, []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if current:
                    seqs[current] = "".join(chunks)
                current, chunks = parse_id(line[1:]), []
            else:
                chunks.append(line)
    if current:
        seqs[current] = "".join(chunks)
    return seqs


def sequence_hash(seq):
    # Same key as biotrainer_core.functions.hashing.calculate_sequence_hash
    return hashlib.sha256(f"{seq}_{len(seq)}".encode()).hexdigest()


def h5_keys(path):
    with h5py.File(path, "r") as f:
        return set(f.keys())


def split_predictor(prefix, h5_base):
    # (id, target, set, sequence) for the extremophile splits; labels and sequences come from the split FASTAs
    entries = []
    for split, set_name in (("train", "train"), ("val", "val"), ("test", "test")):
        with open(SPLITS / f"{prefix}_{split}.faa") as f:
            header, chunks = None, []
            for line in list(f) + [">"]:
                line = line.strip()
                if line.startswith(">"):
                    if header:
                        parts = header.split("||")
                        entries.append((parts[1].split()[0], parts[2].strip(), set_name, "".join(chunks)))
                    header, chunks = line[1:], []
                else:
                    chunks.append(line)
    sources = {e: [f"embeddings/{h5_base}{EMBEDDERS[e]['suffix']}.h5"] for e in EMBEDDERS if EMBEDDERS[e]["suffix"] is not None}
    return entries, sources


def clustered_predictor(split_tsv, binary):
    # Extremophile sequences from the Split 1 FASTAs, non-extremophiles from the filtered Swiss-Prot set; IDs like
    # "Q5HDB2||Q5HDB2" are shortened to "Q5HDB2" in the biotrainer FASTA
    extrem = {}
    for split in ("train", "val", "test"):
        extrem.update(read_faa(SPLITS / f"rep_split1_{split}.faa", lambda h: h.split("||")[1].split()[0]))
    nonextrem = read_faa(NONEXTREM_FAA, lambda h: h.split()[0])
    entries = []
    with open(split_tsv) as f:
        next(f)
        for line in f:
            protein_id, split, label = line.rstrip("\n").split("\t")
            seq = nonextrem[protein_id] if label == "non_extremophile" else extrem[protein_id]
            target = label if not binary or label == "non_extremophile" else "extremophile"
            entries.append((protein_id, target, split, seq))
    sources = {e: [f"embeddings/split1_30pct{EMBEDDERS[e]['suffix']}.h5", f"embeddings/nonextremophile{EMBEDDERS[e]['suffix']}.h5"]
               for e in EMBEDDERS if EMBEDDERS[e]["suffix"] is not None}
    return entries, sources


def write_run(predictor, entries, sources, embedders):
    available = None
    for paths in sources.values():
        keys = set().union(*(h5_keys(p) for p in paths))
        available = keys if available is None else available & keys
    kept = [e for e in entries if e[0] in available]
    print(f"{predictor}: {len(kept):,} of {len(entries):,} proteins have all embeddings")

    run_dir = RUNS_DIR / predictor
    run_dir.mkdir(parents=True, exist_ok=True)
    with open(run_dir / "input.fasta", "w") as f:
        for protein_id, target, set_name, seq in kept:
            f.write(f">{protein_id.split('||')[0]} TARGET={target} SET={set_name}\n{seq}\n")

    EMBED_DIR.mkdir(parents=True, exist_ok=True)
    for embedder in embedders:
        out = EMBED_DIR / f"{predictor}_{embedder}.h5"
        if not out.exists():
            with h5py.File(out, "w") as h5_out:
                handles = [h5py.File(p, "r") for p in sources.get(embedder, [])]
                for protein_id, _, _, seq in kept:
                    key = sequence_hash(seq)
                    if key in h5_out:
                        continue  # identical sequence already stored
                    if embedder == "composition":
                        vector = aa_composition(seq)
                    else:
                        vector = next(h[protein_id][:] for h in handles if protein_id in h)
                    h5_out.create_dataset(key, data=np.asarray(vector, dtype=np.float32))
                    h5_out[key].attrs["original_id"] = protein_id.split("||")[0]
                for h in handles:
                    h.close()
            print(f"  {out}")
        config_dir = run_dir / embedder
        config_dir.mkdir(exist_ok=True)
        (config_dir / "config.yml").write_text(
            CONFIG.format(predictor=predictor, embedder=embedder, embedder_name=EMBEDDERS[embedder]["name"]))


if __name__ == "__main__":
    for predictor, (prefix, h5_base) in SPLIT_PREDICTORS.items():
        write_run(predictor, *split_predictor(prefix, h5_base), list(EMBEDDERS))
    write_run("binary_clustered", *clustered_predictor(SPLITS / "clustered_binary/binary_split.tsv", True), list(EMBEDDERS))
    # The six-class model was only trained with ProtT5
    write_run("sixclass_clustered", *clustered_predictor(SPLITS / "clustered_binary/sixclass_split.tsv", False), ["prott5"])
    print(f"\nConfigs written to {RUNS_DIR}/<predictor>/<embedder>/config.yml")
