# Homology baseline on the algae: share of filtered algae sequences with a prokaryote extremophile hit at different
# E-value thresholds. Identical sequences are counted once.

import csv
import json
from collections import defaultdict
from pathlib import Path

BASE_DIR   = Path(__file__).resolve().parents[2]
HITS_TSV   = BASE_DIR / "results" / "homology_inference" / "homology_hits.tsv"
QUERY_FAA  = BASE_DIR / "results" / "homology_inference" / "filtered_eukaryote_query.faa"
LABELS_TSV = BASE_DIR / "embeddings" / "eukaryote_algae_labels.tsv"
BLACKLIST  = BASE_DIR / "results" / "eukaryote_leakage" / "flagged_protein_ids.txt"
OUT_JSON   = BASE_DIR / "results" / "homology_inference" / "homology_inference_result.json"

THRESHOLDS = [1e-1, 1e-3, 1e-5, 1e-10, 1e-20, 1e-50, 1e-100]


def load_blacklist():
    with open(BLACKLIST) as f:
        return {line.strip() for line in f if line.strip()}


def load_labels(blacklist):
    labels = {}
    with open(LABELS_TSV) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            pid = row["protein_id"]
            if pid not in blacklist:
                labels[pid] = row["broad_class"]
    return labels


def read_sequences():
    seqs, current = {}, None
    with open(QUERY_FAA) as f:
        for line in f:
            if line.startswith(">"):
                current = line[1:].split()[0]
                seqs[current] = []
            else:
                seqs[current].append(line.strip())
    return {protein_id: "".join(chunks) for protein_id, chunks in seqs.items()}


def unique_sequences(labels):
    # One protein per distinct sequence; identical sequences have identical hits
    seen, kept = set(), {}
    for protein_id, seq in read_sequences().items():
        if protein_id in labels and seq not in seen:
            seen.add(seq)
            kept[protein_id] = labels[protein_id]
    return kept


def load_best_evalue(hits_tsv):
    best = {}
    with open(hits_tsv) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue
            query_id, _, _, evalue_str, *_ = parts
            try:
                evalue = float(evalue_str)
            except ValueError:
                continue
            if query_id not in best or evalue < best[query_id]:
                best[query_id] = evalue
    return best


def evaluate(labels, best_evalue, threshold):
    classes = sorted(set(labels.values()))
    total   = defaultdict(int)
    hits    = defaultdict(int)

    for pid, cls in labels.items():
        total[cls] += 1
        if best_evalue.get(pid, float("inf")) <= threshold:
            hits[cls] += 1

    result = {}
    for cls in classes:
        n  = total[cls]
        h  = hits[cls]
        result[cls] = {"recall": h / n if n > 0 else 0.0, "n": n, "hits": h}

    n_total = sum(total.values())
    n_hits  = sum(hits.values())
    result["overall"] = {
        "recall": n_hits / n_total if n_total > 0 else 0.0,
        "n":      n_total,
        "hits":   n_hits,
    }
    return result


if __name__ == "__main__":
    blacklist   = load_blacklist()
    labels      = unique_sequences(load_labels(blacklist))
    best_evalue = load_best_evalue(HITS_TSV)

    print(f"{len(labels):,} unique filtered sequences  |  {len(best_evalue):,} proteins with at least one hit\n")

    results_by_threshold = {}
    for thr in THRESHOLDS:
        res = evaluate(labels, best_evalue, thr)
        results_by_threshold[str(thr)] = res
        r = res["overall"]
        print(f"e-value <= {thr:.0e}:  recall = {r['recall']:.3f}  ({r['hits']:,} / {r['n']:,})")

    out = {
        "thresholds":  [str(t) for t in THRESHOLDS],
        "results":     results_by_threshold,
        "n_sequences": len(labels),
        "n_with_hits": len(best_evalue),
    }

    OUT_JSON.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT_JSON, "w") as f:
        json.dump(out, f, indent=2)
    print(f"\nSaved: {OUT_JSON}")
