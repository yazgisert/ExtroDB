# Runs the pre-trained THERMO model (Koblitz et al. 2025, Pfam-based random forest from bacdive-AI) on our
# proteomes, using pyhmmer for the Pfam annotation. Setup of Pfam-A.hmm and bacdive-AI: see README.md.

import argparse
import csv
import json
import os
import pickle
import sys
import time
from pathlib import Path

import numpy as np
from sklearn.metrics import (
    accuracy_score,
    balanced_accuracy_score,
    f1_score,
    matthews_corrcoef,
    precision_score,
    recall_score,
)

try:
    import pyhmmer
    from pyhmmer.easel import SequenceFile
    from pyhmmer.plan7 import HMMFile
except ImportError:
    sys.exit(
        "pyhmmer is not installed. Run:  pip install pyhmmer\n"
        "Then re-run this script."
    )


EVALUE_THRESHOLD = 1e-20   # same cut-off used in predict.py
N_BOOTSTRAP      = 1000
BOOTSTRAP_SEED   = 42


def load_index(index_path: Path, fasta_base: Path):
    rows = []
    with open(index_path, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            # The algae are evaluated separately with their own thermophile labels (run_thermo_eukaryote_eval.py)
            if row.get("domain_scope") != "prokaryote":
                continue
            # Use the first available FASTA path
            fasta_rel = row.get("local_fasta_path", "").strip()
            if not fasta_rel:
                # fall back to first entry in local_fasta_paths (pipe-separated)
                fasta_rel = row.get("local_fasta_paths", "").split("|")[0].strip()
            if not fasta_rel:
                print(f"  [WARN] No FASTA path for {row['proteome_id']}, skipping.")
                continue

            fasta_abs = fasta_base / fasta_rel
            if not fasta_abs.exists():
                print(f"  [WARN] FASTA not found: {fasta_abs}, skipping.")
                continue

            label = 1 if row.get("thermophile", "no").strip().lower() == "yes" else 0
            rows.append(
                {
                    "proteome_id": row["proteome_id"],
                    "taxon_name":  row.get("taxon_name", ""),
                    "label":       label,
                    "fasta_path":  fasta_abs,
                }
            )
    return rows


def load_thermo_model(model_path: Path):
    with open(model_path, "rb") as f:
        dump = pickle.load(f)
    model      = dump.get("model")
    categories = dump.get("categories")   # list of Pfam accession strings
    if model is None or categories is None:
        sys.exit(f"Unexpected pickle structure in {model_path}. "
                 "Expected keys 'model' and 'categories'.")
    print(f"  Model loaded — {len(categories)} Pfam categories.")
    return model, categories


def load_relevant_hmms(pfam_path: Path, categories: list[str]):
    cat_set = set(categories)
    kept = []
    print(f"  Scanning Pfam-A.hmm for {len(cat_set)} relevant profiles ...", flush=True)
    with HMMFile(str(pfam_path)) as hf:
        for hmm in hf:
            acc  = hmm.accession if isinstance(hmm.accession, str) else (hmm.accession.decode() if hmm.accession else "")
            name = hmm.name if isinstance(hmm.name, str) else (hmm.name.decode() if hmm.name else "")
            # Pfam accessions look like PF00001.23 — strip version suffix
            acc_base = acc.split(".")[0]
            if acc_base in cat_set or acc in cat_set or name in cat_set:
                kept.append(hmm)
    print(f"  Retained {len(kept)} / {len(cat_set)} profiles.")
    return kept


def scan_proteome(fasta_path: Path, hmms: list, evalue: float) -> set:
    hits_pfams = set()
    with SequenceFile(str(fasta_path), digital=True) as sf:
        sequences = sf.read_block()

    for top_hits in pyhmmer.hmmscan(sequences, hmms, cpus=1):
        for hit in top_hits:
            if hit.evalue > evalue:
                continue
            raw = hit.accession if hit.accession else hit.name
            acc = raw if isinstance(raw, str) else raw.decode()
            hits_pfams.add(acc.split(".")[0])   # strip version
    return hits_pfams


def build_feature_vector(pfam_hits: set, categories: list[str]) -> list[int]:
    return [1 if cat in pfam_hits else 0 for cat in categories]


def compute_metrics(y_true, y_pred):
    return {
        "accuracy":          float(accuracy_score(y_true, y_pred)),
        "balanced-accuracy": float(balanced_accuracy_score(y_true, y_pred)),
        "precision":         float(precision_score(y_true, y_pred, zero_division=0)),
        "recall":            float(recall_score(y_true, y_pred, zero_division=0)),
        "f1_score":          float(f1_score(y_true, y_pred, zero_division=0)),
        "macro-f1_score":    float(f1_score(y_true, y_pred, average="macro", zero_division=0)),
        "matthews-corr-coeff": float(matthews_corrcoef(y_true, y_pred)),
    }


def bootstrap_ci(y_true, y_pred, n=N_BOOTSTRAP, seed=BOOTSTRAP_SEED):
    rng = np.random.default_rng(seed)
    n_samples = len(y_true)
    y_true = np.array(y_true)
    y_pred = np.array(y_pred)

    records = []
    for _ in range(n):
        idx = rng.integers(0, n_samples, size=n_samples)
        records.append(compute_metrics(y_true[idx], y_pred[idx]))

    keys = list(records[0].keys())
    result = {}
    for k in keys:
        vals = sorted(r[k] for r in records)
        lower_i = int(0.025 * n)
        upper_i = int(0.975 * n)
        result[k] = {
            "mean":  float(np.mean(vals)),
            "lower": float(vals[lower_i]),
            "upper": float(vals[upper_i]),
        }
    return result


def per_class_metrics(y_true, y_pred, classes):
    out = {}
    y_true = np.array(y_true)
    y_pred = np.array(y_pred)
    for cls_name, cls_val in classes.items():
        mask_true = (y_true == cls_val)
        tp = int(np.sum((y_pred == cls_val) & mask_true))
        fp = int(np.sum((y_pred == cls_val) & ~mask_true))
        fn = int(np.sum((y_pred != cls_val) & mask_true))
        prec = tp / (tp + fp) if (tp + fp) > 0 else 0.0
        rec  = tp / (tp + fn) if (tp + fn) > 0 else 0.0
        f1   = 2 * prec * rec / (prec + rec) if (prec + rec) > 0 else 0.0
        out[cls_name] = {
            "precision": round(prec, 4),
            "recall":    round(rec, 4),
            "f1_score":  round(f1, 4),
            "support":   int(np.sum(mask_true)),
        }
    return out


def main():
    parser = argparse.ArgumentParser(description="Apply THERMO model to extremophile proteomes.")
    parser.add_argument("--pfam",      required=True,  help="Path to Pfam-A.hmm (uncompressed)")
    parser.add_argument("--model",     required=True,  help="Path to thermophile_data.p (pickle)")
    parser.add_argument("--index",     default="searchable_resource/extremophile_master_index.tsv",
                        help="Master index TSV")
    parser.add_argument("--fasta_dir", default=".",    help="Base directory for relative FASTA paths")
    parser.add_argument("--out",       default="results/training_results/thermo_bacdive_result.json",
                        help="Output JSON path")
    parser.add_argument("--evalue",    type=float, default=EVALUE_THRESHOLD,
                        help=f"E-value threshold (default {EVALUE_THRESHOLD})")
    args = parser.parse_args()

    pfam_path   = Path(args.pfam)
    model_path  = Path(args.model)
    index_path  = Path(args.index)
    fasta_base  = Path(args.fasta_dir)
    out_path    = Path(args.out)

    for p, label in [(pfam_path, "Pfam-A.hmm"), (model_path, "model pickle"),
                     (index_path, "master index")]:
        if not p.exists():
            sys.exit(f"ERROR: {label} not found: {p}")

    out_path.parent.mkdir(parents=True, exist_ok=True)

    print("\nTHERMO model evaluation\n")

    print("Step 1/5  Loading master index ...")
    entries = load_index(index_path, fasta_base)
    n_thermo  = sum(e["label"] == 1 for e in entries)
    n_nonther = sum(e["label"] == 0 for e in entries)
    print(f"  {len(entries)} proteomes ({n_thermo} thermophile, {n_nonther} non-thermophile)")

    print("\nStep 2/5  Loading THERMO pickle ...")
    model, categories = load_thermo_model(model_path)

    print("\nStep 3/5  Loading Pfam HMMs ...")
    hmms = load_relevant_hmms(pfam_path, categories)

    print(f"\nStep 4/5  Scanning {len(entries)} proteomes with pyhmmer ...")
    results_per_proteome = []
    t0 = time.time()
    for i, entry in enumerate(entries, 1):
        pfam_hits = scan_proteome(entry["fasta_path"], hmms, args.evalue)
        feat_vec  = build_feature_vector(pfam_hits, categories)
        pred      = int(model.predict([feat_vec])[0])
        try:
            proba = float(model.predict_proba([feat_vec])[0][1])
        except AttributeError:
            proba = float(pred)

        elapsed = time.time() - t0
        eta = (elapsed / i) * (len(entries) - i)
        print(
            f"  [{i:3d}/{len(entries)}] {entry['proteome_id']:<30s} "
            f"label={entry['label']} pred={pred} proba={proba:.3f}  "
            f"pfam_hits={len(pfam_hits)}  ETA {eta:.0f}s",
            flush=True,
        )
        results_per_proteome.append(
            {
                "proteome_id": entry["proteome_id"],
                "taxon_name":  entry["taxon_name"],
                "true_label":  entry["label"],
                "predicted":   pred,
                "proba_thermophile": proba,
                "n_pfam_hits": len(pfam_hits),
            }
        )

    print("\nStep 5/5  Computing metrics ...")
    y_true = [r["true_label"]  for r in results_per_proteome]
    y_pred = [r["predicted"]   for r in results_per_proteome]

    metrics     = compute_metrics(y_true, y_pred)
    metrics_ci  = bootstrap_ci(y_true, y_pred)
    per_class   = per_class_metrics(
        y_true, y_pred,
        classes={"thermophile": 1, "non_thermophile": 0},
    )

    output = {
        "experiment":    "thermo_bacdive",
        "model_source":  "Koblitz et al. 2025 (bacdive-AI, thermophile_data.p)",
        "feature_type":  "Pfam_presence_absence",
        "evalue_threshold": args.evalue,
        "n_proteomes":   len(entries),
        "n_thermophile": n_thermo,
        "n_non_thermophile": n_nonther,
        "n_pfam_categories": len(categories),
        "n_pfam_hmms_matched": len(hmms),
        "note": (
            "All 238 proteomes evaluated with pre-trained external model. "
            "No train/test split — this is a transfer/generalisation test. "
            "Results are NOT directly comparable to published THERMO numbers "
            "(different test set and class distribution)."
        ),
        "test_metrics": {k: round(v, 4) for k, v in metrics.items()},
        "test_metrics_bootstrap": {
            k: {kk: round(vv, 4) for kk, vv in v.items()}
            for k, v in metrics_ci.items()
        },
        "per_class_metrics": per_class,
        "per_proteome_predictions": results_per_proteome,
    }

    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(output, f, indent=2)

    print(f"\nResults saved to {out_path}")
    print("\nSummary:")
    print(f"  Accuracy:   {metrics['accuracy']:.4f}")
    print(f"  Macro-F1:   {metrics['macro-f1_score']:.4f}")
    print(f"  F1 (thermo):{per_class['thermophile']['f1_score']:.4f}")
    print(f"  MCC:        {metrics['matthews-corr-coeff']:.4f}")


if __name__ == "__main__":
    main()
