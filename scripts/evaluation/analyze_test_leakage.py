# Leakage check from the test-vs-train hits: highest identity per test protein, binned, and the IDs with >= 30%
# identity (--clustered for the clustered split).

import argparse
from collections import defaultdict
from pathlib import Path


IDENTITY_THRESHOLD = 30.0
BINS = [(30, 50), (50, 70), (70, 90), (90, 101)]


def group(label):
    return "non_extremophile" if label == "non_extremophile" else "extremophile"


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--clustered", action="store_true", help="analyse the search for the clustered binary split")
    OUT_DIR = Path("results/test_vs_train_search_clustered" if parser.parse_args().clustered else "results/test_vs_train_search")
    HITS    = OUT_DIR / "hits.tsv"
    LABELS  = OUT_DIR / "labels.tsv"
    FLAGGED = OUT_DIR / "flagged_test_ids_30pct.txt"
    SUMMARY = OUT_DIR / "leakage_summary.txt"

    test_group = {}
    with open(LABELS) as f:
        next(f)
        for line in f:
            protein_id, which, label = line.rstrip("\n").split("\t")
            if which == "test":
                test_group[protein_id] = group(label)

    best_identity = defaultdict(float)
    with open(HITS) as f:
        for line in f:
            query, _, pident = line.split("\t")[:3]
            best_identity[query] = max(best_identity[query], float(pident))

    lines = [f"Test proteins vs. training proteins (MMseqs2), highest identity per test protein\n"]
    flagged = []
    for name in ("extremophile", "non_extremophile"):
        ids = [i for i, g in test_group.items() if g == name]
        over = [i for i in ids if best_identity.get(i, 0) >= IDENTITY_THRESHOLD]
        flagged += over
        lines.append(f"{name}: {len(over):,} / {len(ids):,} ({len(over) / len(ids):.1%}) with a hit >= {IDENTITY_THRESHOLD:.0f}% identity")
        for low, high in BINS:
            n = sum(low <= best_identity.get(i, 0) < high for i in ids)
            lines.append(f"    {low:>3}-{min(high, 100):<3}%  {n:>7,}")
        lines.append(f"    no hit   {sum(i not in best_identity for i in ids):>7,}")

    text = "\n".join(lines)
    print(text)
    SUMMARY.write_text(text + "\n")
    FLAGGED.write_text("\n".join(flagged) + "\n")
    print(f"\nSummary saved to {SUMMARY}")
    print(f"Flagged IDs saved to {FLAGGED}")
