"""
write_split_fastas.py

Reads the combined prokaryote FASTA and the stratified split assignments,
then writes six output FASTA files:

  split1_train.faa   split1_val.faa   split1_test.faa
  split2_train.faa   split2_val.faa   split2_test.faa

Protein headers are enriched with the phenotype label:
  >{proteome_id}||{original_protein_id}||{label}

Input files (same folder as this script)
-----------------------------------------
  split_input_prokaryote_singlelabel.faa
  split1_stratified_assignments.tsv
  split2_stratified_assignments.tsv
"""

from __future__ import annotations

import csv
from pathlib import Path

# ── Paths ────────────────────────────────────────────────────────────────────

SPLITS_DIR   = Path(__file__).resolve().parent
COMBINED_FAA = SPLITS_DIR / "split_input_prokaryote_singlelabel.faa"
S1_ASSIGN    = SPLITS_DIR / "split1_stratified_assignments.tsv"
S2_ASSIGN    = SPLITS_DIR / "split2_stratified_assignments.tsv"

PARTITIONS   = ["train", "val", "test"]

# ── Helpers ──────────────────────────────────────────────────────────────────

def read_assignments(path: Path) -> dict[str, dict]:
    """
    Returns dict: proteome_id -> {"label": str, "partition": str}
    """
    result = {}
    with path.open(encoding="utf-8", newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            result[row["proteome_id"]] = {
                "label":     row["label"],
                "partition": row["partition"],
            }
    return result


# ── Main ─────────────────────────────────────────────────────────────────────

def main() -> None:
    print("Reading split assignments ...")
    s1 = read_assignments(S1_ASSIGN)
    s2 = read_assignments(S2_ASSIGN)
    print(f"  Split 1: {len(s1)} organisms")
    print(f"  Split 2: {len(s2)} organisms")

    # Open all 6 output files at once
    handles: dict[tuple[str, str], object] = {}
    for split_name, assign in [("split1", s1), ("split2", s2)]:
        for partition in PARTITIONS:
            out_path = SPLITS_DIR / f"{split_name}_{partition}.faa"
            handles[(split_name, partition)] = out_path.open("w", encoding="utf-8")
            print(f"  Opened: {out_path.name}")

    print(f"\nStreaming combined FASTA: {COMBINED_FAA.name}")
    proteins_written = {("split1", p): 0 for p in PARTITIONS}
    proteins_written.update({("split2", p): 0 for p in PARTITIONS})

    current_proteome = None
    s1_info = None
    s2_info = None
    seq_lines: list[str] = []
    current_header_s1 = None
    current_header_s2 = None

    def flush(s1_header, s2_header, seq_lines, s1_info, s2_info):
        if not seq_lines or s1_info is None:
            return
        seq_block = "".join(seq_lines)
        # Write to Split 1 file
        p1 = s1_info["partition"]
        handles[("split1", p1)].write(s1_header + "\n" + seq_block)
        proteins_written[("split1", p1)] += 1
        # Write to Split 2 file
        p2 = s2_info["partition"]
        handles[("split2", p2)].write(s2_header + "\n" + seq_block)
        proteins_written[("split2", p2)] += 1

    with COMBINED_FAA.open(encoding="utf-8") as fh:
        for i, line in enumerate(fh):
            if i % 500_000 == 0 and i > 0:
                print(f"  {i:,} lines processed ...")
            line = line.rstrip("\n")
            if line.startswith(">"):
                # Flush previous protein
                flush(current_header_s1, current_header_s2, seq_lines,
                      s1_info, s2_info)
                seq_lines = []

                # Parse: >proteome_id||original_protein_id
                raw = line[1:]
                proteome_id = raw.split("||")[0]

                if proteome_id != current_proteome:
                    current_proteome = proteome_id
                    s1_info = s1.get(proteome_id)
                    s2_info = s2.get(proteome_id)

                if s1_info is None:
                    # Organism not in assignments (shouldn't happen)
                    current_header_s1 = None
                    current_header_s2 = None
                else:
                    current_header_s1 = f">{raw}||{s1_info['label']}"
                    current_header_s2 = f">{raw}||{s2_info['label']}"
            else:
                seq_lines.append(line + "\n")

    # Flush last protein
    flush(current_header_s1, current_header_s2, seq_lines, s1_info, s2_info)

    # Close all file handles
    for fh in handles.values():
        fh.close()

    print("\n── Proteins written per file ──")
    total = 0
    for split_name in ["split1", "split2"]:
        for partition in PARTITIONS:
            n = proteins_written[(split_name, partition)]
            fname = f"{split_name}_{partition}.faa"
            print(f"  {fname:<30} {n:>7,} proteins")
            total += n
    print(f"\n  Total writes: {total:,}  "
          f"(= 2 × 667,659 = {2*667659:,} expected)")

    print("\nDone.")


if __name__ == "__main__":
    main()
