# Removes non-extremophile proteins that are near-identical (>= 90% identity, MMseqs2) to any extremophile
# representative and writes nonextremophile_data/nonextremophile_filtered.faa.

from __future__ import annotations

import csv
import subprocess
import tempfile
from pathlib import Path


# 90%, not 30%: homologs across the two classes are expected (ribosomes, polymerases, ...) and are what the
# model has to tell apart; only near-identical copies would be label noise.
IDENTITY_THRESHOLD = 0.90
SENSITIVITY        = 7.5

BASE_DIR   = Path(__file__).resolve().parents[2]
SPLITS_DIR = BASE_DIR / "dataset_splits"

QUERY_FAA   = BASE_DIR / "nonextremophile_data" / "nonextremophile.faa"
OUT_FAA     = BASE_DIR / "nonextremophile_data" / "nonextremophile_filtered.faa"
OUT_HITS    = BASE_DIR / "nonextremophile_data" / "leakage_hits.tsv"
OUT_SUMMARY = BASE_DIR / "nonextremophile_data" / "leakage_summary.txt"

# All extremophile representative FAAs (split1_30pct covers all classes)
EXTREMOPHILE_FAAS = [
    SPLITS_DIR / "rep_split1_train.faa",
    SPLITS_DIR / "rep_split1_val.faa",
    SPLITS_DIR / "rep_split1_test.faa",
]


def combine_fastas(paths: list[Path], out_path: Path) -> int:
    count = 0
    with out_path.open("w") as out:
        for p in paths:
            with p.open(errors="replace") as fh:
                for line in fh:
                    out.write(line)
                    if line.startswith(">"):
                        count += 1
    return count


def count_sequences(faa_path: Path) -> int:
    with open(faa_path) as f:
        return sum(1 for l in f if l.startswith(">"))


def run_mmseqs_search(query: Path, ref: Path, results: Path, tmp_dir: Path) -> None:
    cmd = [
        "mmseqs", "easy-search",
        str(query), str(ref), str(results), str(tmp_dir),
        "--min-seq-id", str(IDENTITY_THRESHOLD),
        "-s", str(SENSITIVITY),
        "--format-output", "query,target,pident,evalue,bits",
        "--threads", "4",
    ]
    print(f"  Running: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(result.stderr[-2000:])
        raise RuntimeError("MMseqs2 failed — see stderr above.")
    for line in result.stderr.splitlines():
        if line.strip():
            print(f"    [mmseqs] {line}")


def parse_query_faa(faa_path: Path) -> dict[str, str]:
    records = {}
    current_id    = None
    current_lines = []

    with open(faa_path) as f:
        for line in f:
            if line.startswith(">"):
                if current_id:
                    records[current_id] = "".join(current_lines)
                # Use first whitespace token — same as what MMseqs2 uses
                current_id    = line[1:].split()[0]
                current_lines = [line]
            else:
                current_lines.append(line)
    if current_id:
        records[current_id] = "".join(current_lines)

    return records


def main() -> None:
    if not QUERY_FAA.exists():
        raise SystemExit(f"Query file not found: {QUERY_FAA}\nRun download_nonextremophile_uniprot.py first.")

    for f in EXTREMOPHILE_FAAS:
        if not f.exists():
            raise SystemExit(f"Missing extremophile FAA: {f}")

    n_query = count_sequences(QUERY_FAA)
    print(f"Non-extremophile sequences to check: {n_query:,}")

    with tempfile.TemporaryDirectory() as tmp_str:
        tmp_dir = Path(tmp_str)
        ref_faa  = tmp_dir / "extremophile_ref.faa"
        raw_hits = tmp_dir / "mmseqs_hits.tsv"

        print("\nCombining extremophile reference FAAs ...")
        n_ref = combine_fastas(EXTREMOPHILE_FAAS, ref_faa)
        print(f"  {n_ref:,} extremophile sequences (reference)")

        print(f"\nRunning MMseqs2 (flagging hits ≥ {IDENTITY_THRESHOLD*100:.0f}% identity) ...")
        mmseqs_tmp = tmp_dir / "mmseqs_tmp"
        mmseqs_tmp.mkdir()
        run_mmseqs_search(QUERY_FAA, ref_faa, raw_hits, mmseqs_tmp)

        flagged_ids: set[str] = set()
        hits = []
        with raw_hits.open() as fh:
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 5:
                    continue
                query_id, target, pident, evalue, bits = parts[:5]
                flagged_ids.add(query_id)
                hits.append((query_id, target, float(pident), evalue, bits))

        print(f"\n  {len(flagged_ids):,} non-extremophile sequences flagged as leakage ({len(flagged_ids)/n_query*100:.2f}%)")

        OUT_HITS.parent.mkdir(exist_ok=True)
        with OUT_HITS.open("w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t")
            w.writerow(["non_extremophile_protein", "extremophile_protein", "pident_%", "evalue", "bits"])
            for row in sorted(hits, key=lambda x: -x[2]):
                w.writerow([row[0], row[1], f"{row[2]:.1f}", row[3], row[4]])
        print(f"  Hit table written: {OUT_HITS}")

        print("\nWriting filtered non-extremophile FAA ...")
        records = parse_query_faa(QUERY_FAA)
        kept = {k: v for k, v in records.items() if k not in flagged_ids}

        with OUT_FAA.open("w") as f:
            for block in kept.values():
                f.write(block)

        n_removed = n_query - len(kept)
        print(f"  Removed:  {n_removed:,} sequences")
        print(f"  Kept:     {len(kept):,} sequences")
        print(f"  Saved to: {OUT_FAA}")

    with OUT_SUMMARY.open("w") as f:
        f.write(f"Non-extremophile leakage check summary\n")
        f.write(f"======================================\n")
        f.write(f"Identity threshold:         {IDENTITY_THRESHOLD*100:.0f}%\n")
        f.write(f"Query sequences:            {n_query:,}\n")
        f.write(f"Reference extremophile seqs:{n_ref:,}\n")
        f.write(f"Flagged and removed:        {len(flagged_ids):,} ({len(flagged_ids)/n_query*100:.2f}%)\n")
        f.write(f"Clean sequences remaining:  {len(kept):,}\n")
        f.write(f"Output:                     {OUT_FAA}\n")

    print("\nSummary:")
    print(f"  Query:              {n_query:,} sequences")
    print(f"  Flagged/removed:    {len(flagged_ids):,} ({len(flagged_ids)/n_query*100:.2f}%)")
    print(f"  Clean remaining:    {len(kept):,}")


if __name__ == "__main__":
    main()
