# Packages every biotrainer predictor as one zip for Zenodo: the input FASTA, and per embedder a ready-to-run config.yml
# plus biotrainer's output (out.yml, log, model checkpoint). Embeddings are not included; the configs name the embedder
# instead, so biotrainer computes them itself. Local absolute paths are removed from the copied files. The algae test
# set is converted to the same FASTA format (all SET=test).
# Run from the repository root: python scripts/biotrainer/package_biotrainer_runs.py

import gzip
import re
import zipfile
from pathlib import Path

RUNS_DIR  = Path("biotrainer_runs")
OUT_DIR   = Path("predictor_resources/biotrainer")
ALGAE_OLD = Path("predictor_resources/zenodo/algae_filtered_test.fasta.gz")
REPO      = str(Path.cwd().resolve())

EMBEDDER_NAMES = {
    "prott5":      "Rostlab/prot_t5_xl_uniref50",
    "esm2_650m":   "facebook/esm2_t33_650M_UR50D",
    "esm2_8m":     "facebook/esm2_t6_8M_UR50D",
    # Mean of biotrainer's one-hot encoding, close to the 20 amino acid frequencies we used
    "composition": "one_hot_encoding",
}

DESCRIPTIONS = {
    "fiveclass_split1_30pct":   "5 extremophile classes, MMseqs2 cluster representatives at 30% identity",
    "fiveclass_split1_20pct":   "5 extremophile classes, MMseqs2 cluster representatives at 20% identity",
    "twelveclass_split2_30pct": "12 extremophile subtypes, MMseqs2 cluster representatives at 30% identity",
    "twelveclass_split2_20pct": "12 extremophile subtypes, MMseqs2 cluster representatives at 20% identity",
    "binary_clustered":         "extremophile vs. non-extremophile; 12,813 non-extremophile cluster representatives (30%) "
                                "and as many extremophiles, equal per class",
    "sixclass_clustered":       "5 extremophile classes plus non-extremophile; test set identical to binary_clustered",
}

README = """{predictor}: {description}

Files
- input.fasta: all proteins in biotrainer format (>ID TARGET=<class> SET=<train|val|test>)
- <embedder>/config.yml: biotrainer config (FNN, one hidden layer of 32, dropout 0.25, Adam, learning rate 1e-3,
  batch size 128, at most 200 epochs, early stopping with patience 10, no class weights, seed 42)
- <embedder>/output/: biotrainer output of our run (out.yml with all test metrics, training log, model checkpoint)

Embedders: {embedders}

Retrain with biotrainer (https://github.com/biocentral/biotrainer, v2.1.1):
    biotrainer train --config prott5/config.yml
biotrainer computes the embeddings itself. Our runs used precomputed embeddings from the biocentral API, which
are not included here because of their size.
"""


def clean(text):
    return text.replace(REPO + "/", "").replace(REPO, ".")


def public_config(config_text, embedder):
    lines = [l for l in config_text.splitlines() if not l.startswith("embeddings_file:") and not l.startswith("# Embeddings")]
    return f"embedder_name: {EMBEDDER_NAMES[embedder]}\n" + "\n".join(lines) + "\n"


def package(predictor_dir):
    predictor = predictor_dir.name
    embedders = sorted(p.parent.name for p in predictor_dir.glob("*/config.yml"))
    zip_path = OUT_DIR / f"{predictor}.zip"
    with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as zf:
        zf.write(predictor_dir / "input.fasta", f"{predictor}/input.fasta")
        zf.writestr(f"{predictor}/README.txt", README.format(predictor=predictor, description=DESCRIPTIONS[predictor],
                                                           embedders=", ".join(embedders)))
        for embedder in embedders:
            run_dir = predictor_dir / embedder
            zf.writestr(f"{predictor}/{embedder}/config.yml", public_config((run_dir / "config.yml").read_text(), embedder))
            outputs = sorted((run_dir / "output").glob("*/out.yml"), key=lambda p: p.stat().st_mtime)
            if not outputs:
                raise SystemExit(f"{run_dir} has no biotrainer output yet")
            for file in sorted(outputs[-1].parent.iterdir()):
                arcname = f"{predictor}/{embedder}/output/{file.name}"
                if file.suffix in (".yml", ".log"):
                    zf.writestr(arcname, clean(file.read_text()))
                else:
                    zf.write(file, arcname)
    print(f"{zip_path}  ({zip_path.stat().st_size / 1e6:.1f} MB, embedders: {', '.join(embedders)})")


def algae_fasta():
    out = OUT_DIR / "algae_filtered_test.fasta.gz"
    with gzip.open(ALGAE_OLD, "rt") as src, gzip.open(out, "wt") as dst:
        for line in src:
            if line.startswith(">"):
                protein_id = line[1:].split()[0]
                target = re.search(r"label=(\S+)", line).group(1)
                dst.write(f">{protein_id} TARGET={target} SET=test\n")
            else:
                dst.write(line)
    print(f"{out}  ({out.stat().st_size / 1e6:.1f} MB)")


if __name__ == "__main__":
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    for predictor_dir in sorted(p for p in RUNS_DIR.iterdir() if (p / "input.fasta").exists()):
        package(predictor_dir)
    algae_fasta()
