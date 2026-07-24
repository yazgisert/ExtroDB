"""
check_algae_leakage.py

Checks sequence similarity between the 13 eukaryotic algae proteomes (held-out
test set) and the balanced prokaryote train + val representatives.

Steps
-----
1. Read the master index to locate all downloaded algae FASTAs.
2. Combine them into a single algae query FASTA.
3. Combine rep_split1_train.faa + rep_split1_val.faa into a prokaryote
   reference FASTA (these are the sequences the model will actually train on).
4. Run MMseqs2 easy-search: algae vs. prokaryote train+val.
5. Report algae proteins with any hit at or above IDENTITY_THRESHOLD.

Output
------
algae_leakage_hits.tsv    — every hit above threshold with identity and e-value
algae_leakage_summary.tsv — per-organism count of flagged proteins

Reproducibility
---------------
IDENTITY_THRESHOLD = 0.30  (matches the MMseqs2 clustering threshold used for
                             the prokaryote split)
RANDOM_SEED = 42
"""

from __future__ import annotations

import csv
import subprocess
import tempfile
from collections import defaultdict
from pathlib import Path

# ── Config ────────────────────────────────────────────────────────────────────

IDENTITY_THRESHOLD = 0.30   # flag algae proteins with ≥30% identity to train/val
SENSITIVITY        = 7.5    # MMseqs2 -s flag; higher = more sensitive / slower

BASE_DIR   = Path(__file__).resolve().parent.parent   # thesis_data_collection/
SPLITS_DIR = Path(__file__).resolve().parent          # dataset_splits/
INDEX_PATH = BASE_DIR / "searchable_resource" / "extremophile_master_index.tsv"

TRAIN_FAA  = SPLITS_DIR / "rep_split1_train.faa"
VAL_FAA    = SPLITS_DIR / "rep_split1_val.faa"

OUT_HITS    = SPLITS_DIR / "algae_leakage_hits.tsv"
OUT_SUMMARY = SPLITS_DIR / "algae_leakage_summary.tsv"

# ── Helpers ───────────────────────────────────────────────────────────────────

def find_algae_fastas() -> list[tuple[str, Path]]:
    """
    Read the master index and return (taxon_name, fasta_path) for every
    downloaded eukaryotic algae organism.
    """
    entries = []
    with INDEX_PATH.open(encoding="utf-8", newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if "eukaryotic_algae" not in row.get("dataset_branch", ""):
                continue
            if row.get("proteome_status") != "downloaded":
                continue
            rel_path = row.get("local_fasta_path", "").strip()
            if not rel_path:
                print(f"  WARNING: no local_fasta_path for {row['taxon_name']} — skipped.")
                continue
            fasta_path = BASE_DIR / rel_path
            if not fasta_path.exists():
                print(f"  WARNING: FASTA not found on disk: {fasta_path} — skipped.")
                continue
            entries.append((row["taxon_name"], fasta_path))
    return entries


def combine_fastas(paths: list[Path], out_path: Path) -> int:
    """Concatenate multiple FASTA files into one. Returns total sequence count."""
    count = 0
    with out_path.open("w", encoding="utf-8") as out:
        for p in paths:
            with p.open(encoding="utf-8", errors="replace") as fh:
                for line in fh:
                    out.write(line)
                    if line.startswith(">"):
                        count += 1
    return count


def run_mmseqs_search(query: Path, ref: Path, results: Path, tmp_dir: Path) -> None:
    """Run MMseqs2 easy-search: every query sequence vs. every reference sequence."""
    cmd = [
        "mmseqs", "easy-search",
        str(query), str(ref), str(results), str(tmp_dir),
        "--min-seq-id", str(IDENTITY_THRESHOLD),
        "-s", str(SENSITIVITY),
        "--format-output", "query,target,pident,evalue,bits",
        "--threads", "4",
    ]
    print(f"  Command: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(result.stderr[-2000:])
        raise RuntimeError("MMseqs2 easy-search failed — see stderr above.")
    for line in result.stderr.splitlines():
        if line.strip():
            print(f"    [mmseqs] {line}")


# ── Main ──────────────────────────────────────────────────────────────────────

def main() -> None:
    # 1. Locate algae FASTAs
    print("Finding algae FASTAs from master index...")
    algae = find_algae_fastas()
    print(f"  {len(algae)} downloaded algae organisms found:")
    for name, path in algae:
        print(f"    {name} — {path.name}")

    if not algae:
        raise SystemExit("No algae FASTAs found — cannot proceed.")

    # Check train/val files exist
    for f in (TRAIN_FAA, VAL_FAA):
        if not f.exists():
            raise SystemExit(
                f"Missing: {f.name}\n"
                "Run assign_splits_mmseqs_representatives.py first."
            )

    with tempfile.TemporaryDirectory() as tmp_str:
        tmp_dir   = Path(tmp_str)
        query_faa = tmp_dir / "algae_query.faa"
        ref_faa   = tmp_dir / "prokaryote_train_val.faa"
        raw_hits  = tmp_dir / "mmseqs_hits.tsv"

        # 2. Combine algae proteins into one query FASTA
        print(f"\nCombining algae FASTAs...")
        n_query = combine_fastas([p for _, p in algae], query_faa)
        print(f"  {n_query:,} algae proteins total")

        # 3. Combine prokaryote train + val into reference FASTA
        print(f"\nCombining prokaryote train + val...")
        n_ref = combine_fastas([TRAIN_FAA, VAL_FAA], ref_faa)
        print(f"  {n_ref:,} prokaryote proteins (reference)")

        # 4. Run MMseqs2
        print(f"\nRunning MMseqs2 easy-search "
              f"(reporting hits ≥ {IDENTITY_THRESHOLD*100:.0f}% identity)...")
        mmseqs_tmp = tmp_dir / "mmseqs_tmp"
        mmseqs_tmp.mkdir()
        run_mmseqs_search(query_faa, ref_faa, raw_hits, mmseqs_tmp)

        # 5. Parse hits
        print("\nParsing results...")
        hits: list[dict] = []
        with raw_hits.open(encoding="utf-8") as fh:
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 5:
                    continue
                query, target, pident, evalue, bits = parts[:5]
                hits.append({
                    "query":  query,
                    "target": target,
                    "pident": float(pident),
                    "evalue": evalue,
                    "bits":   bits,
                })
        print(f"  {len(hits):,} hits at or above {IDENTITY_THRESHOLD*100:.0f}% identity")

        # 6. Write full hit table (sorted by identity descending)
        with OUT_HITS.open("w", encoding="utf-8", newline="") as fh:
            w = csv.writer(fh, delimiter="\t")
            w.writerow(["algae_protein", "prokaryote_protein", "pident_%", "evalue", "bits"])
            for h in sorted(hits, key=lambda x: -x["pident"]):
                w.writerow([h["query"], h["target"], f"{h['pident']:.1f}", h["evalue"], h["bits"]])
        print(f"  Written: {OUT_HITS.name}")

        # 7. Per-organism summary
        # Map proteome_id prefix → organism name
        prefix_to_name: dict[str, str] = {}
        total_per_org:  dict[str, int] = defaultdict(int)
        for name, path in algae:
            with path.open(encoding="utf-8", errors="replace") as fh:
                for line in fh:
                    if line.startswith(">"):
                        # Headers are either plain accessions or proteome_id||protein_id
                        proteome_id = line[1:].split("||")[0].split()[0].strip()
                        prefix_to_name[proteome_id] = name
                        total_per_org[name] += 1

        flagged_per_org: dict[str, set[str]] = defaultdict(set)
        for h in hits:
            proteome_id = h["query"].split("||")[0].split()[0].strip()
            org = prefix_to_name.get(proteome_id, proteome_id)
            flagged_per_org[org].add(h["query"])

        grand_total    = sum(total_per_org.values())
        grand_flagged  = sum(len(v) for v in flagged_per_org.values())

        with OUT_SUMMARY.open("w", encoding="utf-8", newline="") as fh:
            w = csv.writer(fh, delimiter="\t")
            w.writerow(["organism", "total_proteins", "flagged_proteins", "flagged_pct"])
            for name, _ in sorted(algae, key=lambda x: x[0]):
                n_total   = total_per_org.get(name, 0)
                n_flagged = len(flagged_per_org.get(name, set()))
                pct       = n_flagged / n_total * 100 if n_total else 0.0
                w.writerow([name, n_total, n_flagged, f"{pct:.1f}%"])
            w.writerow(["TOTAL", grand_total, grand_flagged,
                        f"{grand_flagged/grand_total*100:.1f}%" if grand_total else "0.0%"])
        print(f"  Written: {OUT_SUMMARY.name}")

    # 8. Print summary
    unique_flagged = len(set(h["query"] for h in hits))
    print(f"\n── Leakage Summary ──────────────────────────────────────")
    print(f"  Algae proteins queried :   {n_query:,}")
    print(f"  Prokaryote train+val ref:  {n_ref:,}")
    print(f"  Identity threshold:        {IDENTITY_THRESHOLD*100:.0f}%")
    print(f"  Unique algae proteins")
    print(f"  with ≥1 hit in train/val:  {unique_flagged:,} "
          f"({unique_flagged/n_query*100:.2f}%)")
    print(f"\nDone. Inspect {OUT_SUMMARY.name} for the per-organism breakdown.")


if __name__ == "__main__":
    main()
