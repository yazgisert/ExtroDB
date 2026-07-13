"""
assign_splits.py

Reads the MMseqs2 cluster output and the organism metadata, then assigns each
organism to train / val / test for two separate splits:

  Split 1  —  coarse label   (phenotype_focus: 5 classes)
  Split 2  —  fine-grained   (split2_label: subtypes)

Algorithm
---------
1. Build a Union-Find over organisms: two organisms are in the same connected
   component if any of their proteins ended up in the same MMseqs2 cluster.
   All organisms in one component MUST go to the same partition (train/val/test)
   to prevent sequence leakage.

2. Greedy assignment: sort components largest-first (by protein count), then
   assign each component to whichever partition keeps the per-label protein
   fractions closest to the 70 / 10 / 20 targets.

3. Write:
     split1_assignments.tsv   — one row per organism (Split 1)
     split2_assignments.tsv   — one row per organism (Split 2)
     split1_summary.tsv       — per-label counts and fractions (Split 1)
     split2_summary.tsv       — per-label counts and fractions (Split 2)

Inputs (in the same folder as this script)
-------------------------------------------
  split_input_metadata.tsv
  mmseqs2_cluster_cluster.tsv
"""

from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path

# ── Paths ────────────────────────────────────────────────────────────────────

SPLITS_DIR  = Path(__file__).resolve().parent
META_TSV    = SPLITS_DIR / "split_input_metadata.tsv"
CLUSTER_TSV = SPLITS_DIR / "mmseqs2_cluster_cluster.tsv"

OUT_S1      = SPLITS_DIR / "split1_assignments.tsv"
OUT_S2      = SPLITS_DIR / "split2_assignments.tsv"
OUT_S1_SUM  = SPLITS_DIR / "split1_summary.tsv"
OUT_S2_SUM  = SPLITS_DIR / "split2_summary.tsv"

# ── Targets ──────────────────────────────────────────────────────────────────

PARTITIONS  = ["train", "val", "test"]
TARGET      = {"train": 0.70, "val": 0.10, "test": 0.20}

# ── Union-Find ───────────────────────────────────────────────────────────────

class UnionFind:
    def __init__(self, items):
        self.parent = {x: x for x in items}
        self.rank   = {x: 0 for x in items}

    def find(self, x):
        while self.parent[x] != x:
            self.parent[x] = self.parent[self.parent[x]]   # path compression
            x = self.parent[x]
        return x

    def union(self, a, b):
        ra, rb = self.find(a), self.find(b)
        if ra == rb:
            return
        if self.rank[ra] < self.rank[rb]:
            ra, rb = rb, ra
        self.parent[rb] = ra
        if self.rank[ra] == self.rank[rb]:
            self.rank[ra] += 1

# ── Helpers ──────────────────────────────────────────────────────────────────

def read_metadata(path: Path) -> dict:
    meta = {}
    with path.open(encoding="utf-8", newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            pid = row["proteome_id"]
            count_str = row.get("protein_count", "0").strip()
            meta[pid] = {
                "split1_label": row["phenotype_focus"],
                "split2_label": row["split2_label"],
                "protein_count": int(count_str) if count_str.isdigit() else 0,
            }
    return meta


def build_components(cluster_tsv: Path, organisms: set) -> dict:
    """
    Parse mmseqs2_cluster_cluster.tsv and return connected components.
    Returns: dict mapping component_root -> list of organism ids.
    """
    uf = UnionFind(organisms)

    print(f"  Parsing {cluster_tsv.name} ...")
    with cluster_tsv.open(encoding="utf-8") as fh:
        for i, line in enumerate(fh):
            if i % 200_000 == 0 and i > 0:
                print(f"    {i:,} lines processed")
            # MMseqs2 outputs tab-separated; fall back to space-split
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                parts = line.rstrip("\n").split()
            if len(parts) < 2:
                continue
            rep_org = parts[0].split("||")[0]
            mem_org = parts[1].split("||")[0]
            # Only link organisms that are in our metadata (ignore stray ids)
            if rep_org in organisms and mem_org in organisms:
                uf.union(rep_org, mem_org)
    print(f"    {i+1:,} lines total")

    components: dict[str, list[str]] = defaultdict(list)
    for org in organisms:
        components[uf.find(org)].append(org)
    return dict(components)


def assign_partitions(
    components: dict,
    meta: dict,
    label_key: str,
) -> dict[str, str]:
    """
    Greedy assignment of connected components to train / val / test.
    Returns: dict proteome_id -> partition name.
    """
    # Total protein count per label (across all organisms)
    total_by_label: dict[str, int] = defaultdict(int)
    for info in meta.values():
        total_by_label[info[label_key]] += info["protein_count"]

    # Running protein count assigned to each partition per label
    current: dict[str, dict[str, int]] = {
        p: defaultdict(int) for p in PARTITIONS
    }

    # Sort components largest-first so big constraints are resolved early
    sorted_comps = sorted(
        components.values(),
        key=lambda orgs: sum(meta[o]["protein_count"] for o in orgs),
        reverse=True,
    )

    assignments: dict[str, str] = {}

    for comp_orgs in sorted_comps:
        # Protein count per label in this component
        comp_by_label: dict[str, int] = defaultdict(int)
        for o in comp_orgs:
            comp_by_label[meta[o][label_key]] += meta[o]["protein_count"]

        # Choose the partition that minimises deviation from targets
        best_partition, best_cost = None, float("inf")
        for partition in PARTITIONS:
            cost = 0.0
            for label, comp_count in comp_by_label.items():
                total = total_by_label[label]
                if total == 0:
                    continue
                new_frac = (current[partition][label] + comp_count) / total
                cost += abs(new_frac - TARGET[partition])
            if cost < best_cost:
                best_cost = cost
                best_partition = partition

        for o in comp_orgs:
            assignments[o] = best_partition
            current[best_partition][meta[o][label_key]] += meta[o]["protein_count"]

    return assignments


def write_assignments(
    assignments: dict,
    meta: dict,
    label_key: str,
    out_path: Path,
) -> None:
    with out_path.open("w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["proteome_id", "label", "partition", "protein_count"])
        for org in sorted(assignments):
            w.writerow([
                org,
                meta[org][label_key],
                assignments[org],
                meta[org]["protein_count"],
            ])
    print(f"  Written: {out_path.name}")


def write_summary(
    assignments: dict,
    meta: dict,
    label_key: str,
    out_path: Path,
) -> None:
    counts: dict = defaultdict(lambda: defaultdict(lambda: {"orgs": 0, "proteins": 0}))
    for org, partition in assignments.items():
        label = meta[org][label_key]
        counts[label][partition]["orgs"]     += 1
        counts[label][partition]["proteins"] += meta[org]["protein_count"]

    with out_path.open("w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["label", "partition", "n_organisms", "n_proteins", "pct_proteins"])
        for label in sorted(counts):
            total = sum(counts[label][p]["proteins"] for p in PARTITIONS)
            for partition in PARTITIONS:
                d = counts[label][partition]
                pct = d["proteins"] / total * 100 if total > 0 else 0.0
                w.writerow([label, partition, d["orgs"], d["proteins"], f"{pct:.1f}%"])
    print(f"  Written: {out_path.name}")


def report_cross_phenotype(components: dict, meta: dict, label_key: str, split_name: str) -> None:
    """Warn if any component spans more than one phenotype class (leakage risk)."""
    cross = []
    for orgs in components.values():
        labels = {meta[o][label_key] for o in orgs}
        if len(labels) > 1:
            cross.append((sorted(orgs), sorted(labels)))
    if cross:
        print(f"  WARNING: {len(cross)} component(s) span multiple {split_name} labels:")
        for orgs, labels in cross:
            print(f"    labels={labels}  organisms={orgs}")
    else:
        print(f"  {split_name}: no cross-label components — clean!")


# ── Main ─────────────────────────────────────────────────────────────────────

def main() -> None:
    print(f"Reading metadata: {META_TSV.name}")
    meta = read_metadata(META_TSV)
    print(f"  {len(meta)} organisms loaded")

    print(f"\nBuilding connected components from MMseqs2 clusters ...")
    components = build_components(CLUSTER_TSV, set(meta.keys()))
    n_single = sum(1 for orgs in components.values() if len(orgs) == 1)
    n_multi  = len(components) - n_single
    print(f"  {len(components)} components total  "
          f"({n_single} singletons, {n_multi} multi-organism)")

    print("\n── Checking for cross-label components ──")
    report_cross_phenotype(components, meta, "split1_label", "Split 1 (coarse)")
    report_cross_phenotype(components, meta, "split2_label", "Split 2 (subtypes)")

    print("\n── Split 1: coarse phenotype (5 classes) ──")
    s1 = assign_partitions(components, meta, "split1_label")
    write_assignments(s1, meta, "split1_label", OUT_S1)
    write_summary(s1, meta, "split1_label", OUT_S1_SUM)

    print("\n── Split 2: fine-grained subtypes ──")
    s2 = assign_partitions(components, meta, "split2_label")
    write_assignments(s2, meta, "split2_label", OUT_S2)
    write_summary(s2, meta, "split2_label", OUT_S2_SUM)

    print("\n── Split 1 summary ──")
    with OUT_S1_SUM.open(encoding="utf-8") as fh:
        print(fh.read())

    print("── Split 2 summary ──")
    with OUT_S2_SUM.open(encoding="utf-8") as fh:
        print(fh.read())

    print("Done.")


if __name__ == "__main__":
    main()
