# Runs THERMO on the 13 algae proteomes (one prediction per species) to test whether the
# prokaryote-trained model recognises the 4 thermophilic algae. Setup: see README.md.

import argparse
import csv
import json
import pickle
import sys
import time
import warnings
from pathlib import Path

import numpy as np

warnings.filterwarnings("ignore", category=UserWarning)   # sklearn version mismatch

try:
    import pyhmmer
    from pyhmmer.easel import SequenceFile
    from pyhmmer.plan7 import HMMFile
except ImportError:
    sys.exit("pyhmmer not found. Run: pip install pyhmmer --break-system-packages")


EVALUE_THRESHOLD = 1e-20   # same cut-off as the original THERMO predict.py
INDEX_PATH       = Path("searchable_resource/extremophile_master_index.tsv")
FASTA_BASE       = Path("algae_raw_protein_fastas")


def load_algae_entries(index_path: Path):
    entries = []
    with open(index_path, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            if "eukaryot" not in row.get("dataset_branch", ""):
                continue

            fasta_rel = (row.get("local_fasta_path", "") or
                         row.get("local_fasta_paths", "").split("|")[0]).strip()
            if not fasta_rel:
                continue

            fasta_abs = FASTA_BASE / Path(fasta_rel).name
            if not fasta_abs.exists():
                print(f"  [SKIP] FASTA not found: {fasta_abs.name}")
                continue

            is_thermo = any(
                row.get(col, "no").strip().lower() == "yes"
                for col in ("thermophilic_alga", "high_temperature_alga",
                            "thermoacidophilic_alga")
            )

            if is_thermo:
                broad_class = "thermophile"
            elif row.get("psychrophilic_alga", "no") == "yes":
                broad_class = "psychrophile"
            elif row.get("saline_adapted_alga", "no") == "yes" or row.get("halophilic_alga", "no") == "yes":
                broad_class = "halophile"
            elif row.get("acidophilic_alga", "no") == "yes":
                broad_class = "acidophile"
            else:
                broad_class = "unknown"

            entries.append({
                "proteome_id": row["proteome_id"],
                "species":     row.get("taxon_name", ""),
                "label":       1 if is_thermo else 0,
                "broad_class": broad_class,
                "fasta_path":  fasta_abs,
            })

    return entries


def load_thermo_model(model_path: Path):
    with open(model_path, "rb") as f:
        dump = pickle.load(f)
    model      = dump["model"]
    categories = list(dump["categories"])
    print(f"  Model loaded  — {len(categories)} Pfam categories, "
          f"{model.n_estimators} trees")
    return model, categories


def load_relevant_hmms(pfam_path: Path, categories: list):
    cat_set = set(categories)
    kept = []
    is_subset = "thermo_profiles" in str(pfam_path)
    if is_subset:
        print(f"  Loading pre-extracted HMM subset from {pfam_path.name} ...", flush=True)
    else:
        print(f"  Scanning {pfam_path.name} for {len(cat_set)} profiles ...", flush=True)
    with HMMFile(str(pfam_path)) as hf:
        for hmm in hf:
            acc  = (hmm.accession or b"").decode() if isinstance(hmm.accession, bytes) else (hmm.accession or "")
            name = (hmm.name or b"").decode()      if isinstance(hmm.name, bytes)      else (hmm.name or "")
            acc_base = acc.split(".")[0]
            if is_subset or acc_base in cat_set or acc in cat_set or name in cat_set:
                kept.append(hmm)
    print(f"  Loaded {len(kept)} profiles.")
    return kept


def scan_proteome(fasta_path: Path, hmms: list, evalue: float) -> set:
    with SequenceFile(str(fasta_path), digital=True) as sf:
        sequences = sf.read_block()
    hits = set()
    for top_hits in pyhmmer.hmmscan(sequences, hmms, cpus=4):
        for hit in top_hits:
            if hit.evalue > evalue:
                continue
            raw = hit.accession if hit.accession else hit.name
            acc = raw.decode() if isinstance(raw, bytes) else raw
            hits.add(acc.split(".")[0])
    return hits


def build_feature_vector(pfam_hits: set, categories: list) -> list:
    return [1 if cat in pfam_hits else 0 for cat in categories]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--pfam",  default="thermo_profiles.hmm",
                        help="Path to Pfam HMM file (default: thermo_profiles.hmm subset)")
    parser.add_argument("--model", required=True,  help="Path to thermophile_data.p")
    parser.add_argument("--out",   default="results/training_results/thermo_eukaryote_eval_result.json")
    parser.add_argument("--evalue", type=float, default=EVALUE_THRESHOLD)
    parser.add_argument("--batch", type=int, default=None,
                        help="Process only this batch (1, 2, or 3). Omit for full run.")
    args = parser.parse_args()

    pfam_path  = Path(args.pfam)
    model_path = Path(args.model)
    out_path   = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    for p, label in [(pfam_path, "Pfam-A.hmm"), (model_path, "model pickle"),
                     (INDEX_PATH, "master index")]:
        if not p.exists():
            sys.exit(f"ERROR: {label} not found: {p}")

    print("\nTHERMO on eukaryotic algae\n")

    print("Step 1/4  Loading algal proteome entries ...")
    all_entries = load_algae_entries(INDEX_PATH)

    # Batch split: Polarella (58K proteins) needs its own batch
    polarella_entries = [e for e in all_entries if "Polarella" in e["species"]]
    medium_entries    = [e for e in all_entries if "Polarella" not in e["species"]
                         and e["fasta_path"].stat().st_size > 7_000_000]
    small_entries     = [e for e in all_entries if "Polarella" not in e["species"]
                         and e["fasta_path"].stat().st_size <= 7_000_000]

    if args.batch == 1:
        entries = polarella_entries
        print("  Batch 1: Polarella (large proteome)")
    elif args.batch == 2:
        entries = medium_entries
        print("  Batch 2: medium proteomes")
    elif args.batch == 3:
        entries = small_entries
        print("  Batch 3: small proteomes")
    else:
        entries = all_entries

    n_thermo = sum(e["label"] == 1 for e in entries)
    n_other  = sum(e["label"] == 0 for e in entries)
    print(f"  {len(entries)} proteomes  ({n_thermo} thermophile, {n_other} non-thermophile)")
    for e in entries:
        print(f"    {e['species']:<35s} label={e['label']}  class={e['broad_class']}")

    print("\nStep 2/4  Loading THERMO model ...")
    model, categories = load_thermo_model(model_path)

    print("\nStep 3/4  Loading Pfam HMMs ...")
    hmms = load_relevant_hmms(pfam_path, categories)

    print(f"\nStep 4/4  Scanning {len(entries)} proteomes ...")
    per_proteome = []
    t0 = time.time()
    for i, entry in enumerate(entries, 1):
        pfam_hits = scan_proteome(entry["fasta_path"], hmms, args.evalue)
        feat_vec  = build_feature_vector(pfam_hits, categories)
        pred      = int(model.predict([feat_vec])[0])
        proba     = float(model.predict_proba([feat_vec])[0][1])

        elapsed = time.time() - t0
        eta     = (elapsed / i) * (len(entries) - i)
        print(
            f"  [{i:2d}/{len(entries)}] {entry['species']:<35s} "
            f"true={entry['label']}  pred={pred}  prob={proba:.3f}  "
            f"pfam_hits={len(pfam_hits)}  ETA {eta:.0f}s",
            flush=True,
        )
        per_proteome.append({
            "proteome_id":      entry["proteome_id"],
            "species":          entry["species"],
            "broad_class":      entry["broad_class"],
            "true_label":       entry["label"],
            "predicted":        pred,
            "proba_thermophile": proba,
            "n_pfam_hits":      len(pfam_hits),
            "correct":          entry["label"] == pred,
        })

    thermo_entries = [r for r in per_proteome if r["true_label"] == 1]
    other_entries  = [r for r in per_proteome if r["true_label"] == 0]

    thermo_recall  = (sum(r["predicted"] == 1 for r in thermo_entries) / len(thermo_entries)
                      if thermo_entries else None)
    other_specificity = (sum(r["predicted"] == 0 for r in other_entries) / len(other_entries)
                         if other_entries else None)

    print("\nResults:")
    print(f"  Thermophile recall (n={len(thermo_entries)}):    "
          f"{thermo_recall:.3f}" if thermo_recall is not None else "  No thermophile entries")
    print(f"  Non-thermophile specificity (n={len(other_entries)}): "
          f"{other_specificity:.3f}" if other_specificity is not None else "  No non-thermophile entries")
    print("  Per proteome:")
    for r in per_proteome:
        status = "CORRECT" if r["correct"] else "WRONG  "
        print(f"    {status}  {r['species']:<35s} true={r['true_label']}  "
              f"pred={r['predicted']}  prob={r['proba_thermophile']:.3f}")

    output = {
        "experiment":        "thermo_eukaryote_eval",
        "model_source":      "Koblitz et al. 2025 (bacdive-AI, thermophile_data.p)",
        "feature_type":      "Pfam_presence_absence",
        "evalue_threshold":  args.evalue,
        "n_proteomes":       len(per_proteome),
        "n_thermophile":     n_thermo,
        "n_non_thermophile": n_other,
        "n_pfam_categories": len(categories),
        "n_pfam_hmms_matched": len(hmms),
        "note": (
            "THERMO is a prokaryote-trained thermophile classifier. "
            "This evaluation tests cross-domain transfer to eukaryote algae. "
            "Classification is at proteome level (one prediction per species). "
            "Thermophile label uses thermophilic_alga / high_temperature_alga / "
            "thermoacidophilic_alga columns — not the prokaryote thermophile column."
        ),
        "thermophile_recall":      round(thermo_recall, 4) if thermo_recall is not None else None,
        "non_thermophile_specificity": round(other_specificity, 4) if other_specificity is not None else None,
        "per_proteome": per_proteome,
    }

    # If running in batch mode, save to a batch-specific file and also try to aggregate
    if args.batch is not None:
        batch_path = out_path.parent / f"thermo_eukaryote_batch{args.batch}.json"
        with open(batch_path, "w", encoding="utf-8") as f:
            json.dump({"per_proteome": per_proteome}, f, indent=2)
        print(f"Batch results saved to {batch_path}")

        # Try to aggregate all batches
        all_batch_files = list(out_path.parent.glob("thermo_eukaryote_batch*.json"))
        all_pp = []
        for bf in sorted(all_batch_files):
            with open(bf) as f:
                all_pp.extend(json.load(f)["per_proteome"])
        if len(all_pp) > 0:
            t_entries = [r for r in all_pp if r["true_label"] == 1]
            o_entries  = [r for r in all_pp if r["true_label"] == 0]
            t_recall   = (sum(r["predicted"] == 1 for r in t_entries) / len(t_entries)
                          if t_entries else None)
            o_spec     = (sum(r["predicted"] == 0 for r in o_entries) / len(o_entries)
                          if o_entries else None)
            output["per_proteome"]  = all_pp
            output["n_proteomes"]   = len(all_pp)
            output["n_thermophile"] = len(t_entries)
            output["n_non_thermophile"] = len(o_entries)
            output["thermophile_recall"] = round(t_recall, 4) if t_recall is not None else None
            output["non_thermophile_specificity"] = round(o_spec, 4) if o_spec is not None else None
            print(f"\nAggregated {len(all_pp)} proteomes from {len(all_batch_files)} batch(es)")
            if t_recall is not None:
                print(f"  Thermophile recall: {t_recall:.3f}  ({len(t_entries)} species)")
            if o_spec is not None:
                print(f"  Non-thermophile specificity: {o_spec:.3f}  ({len(o_entries)} species)")

    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(output, f, indent=2)
    print(f"Results saved to {out_path}")


if __name__ == "__main__":
    main()
