# Compares the biotrainer runs with the earlier results of our own training scripts: test accuracy, macro-F1 and MCC
# with 95% bootstrap CIs per predictor and embedder, and whether the two confidence intervals overlap.
# Run from the repository root (biotrainer env): python scripts/biotrainer/compare_with_previous.py

import json
from pathlib import Path

import yaml

RUNS_DIR = Path("biotrainer_runs")
RESULTS  = Path("results/training_results")
OUT_TSV  = Path("results/biotrainer/biotrainer_comparison.tsv")

SUFFIX = {"prott5": "", "esm2_650m": "_esm2", "esm2_8m": "_esm2_8m", "composition": "_onehot"}
SPLIT_RESULT = {"fiveclass_split1_30pct": "split1_30pct", "fiveclass_split1_20pct": "split1_20pct",
                "twelveclass_split2_30pct": "split2_30pct", "twelveclass_split2_20pct": "split2_20pct"}
BINARY_NAME = {"prott5": "prott5", "esm2_650m": "esm2_650m", "esm2_8m": "esm2_8m", "composition": "onehot"}


def latest_out_yml(run_dir):
    outs = sorted((run_dir / "output").glob("*/out.yml"), key=lambda p: p.stat().st_mtime)
    return outs[-1] if outs else None


def biotrainer_metrics(out_yml):
    test = yaml.safe_load(open(out_yml))["test_results"]["test"]
    boot = {m["name"]: m for m in test["bootstrapped_metrics"]}
    point = test["inference_result"]["metrics"]
    # Multi-class runs report macro-f1_score; binary runs only report f1_score
    f1 = "macro-f1_score" if "macro-f1_score" in point else "f1_score"
    names = {"accuracy": "accuracy", "macro_f1": f1, "mcc": "matthews-corr-coeff"}
    return {k: (point[v], boot[v]["lower"], boot[v]["upper"]) for k, v in names.items()}


def as_triple(entry):
    if isinstance(entry, dict):
        return entry.get("value"), entry.get("ci_lower"), entry.get("ci_upper")
    return entry, None, None


def previous_metrics(predictor, embedder):
    if predictor in SPLIT_RESULT:
        path = RESULTS / f"{SPLIT_RESULT[predictor]}{SUFFIX[embedder]}_result.json"
        r = json.load(open(path))
        m, b = r["test_metrics"], r["test_metrics_bootstrap"]
        names = {"accuracy": "accuracy", "macro_f1": "macro-f1_score", "mcc": "matthews-corr-coeff"}
        return {k: (m[v], b[v]["lower"], b[v]["upper"]) for k, v in names.items()}
    if predictor == "binary_clustered":
        m = json.load(open(RESULTS / f"binary_{BINARY_NAME[embedder]}_clustered_result.json"))["test_metrics"]
        return {"accuracy": as_triple(m["accuracy"]), "macro_f1": as_triple(m["macro_f1"]), "mcc": as_triple(m["mcc"])}
    if predictor == "sixclass_clustered":
        m = json.load(open(RESULTS / "sixclass_clustered_result.json"))["sixclass_metrics"]
        return {"accuracy": as_triple(m["accuracy"]), "macro_f1": as_triple(m["macro_f1"]), "mcc": as_triple(m["mcc"])}
    return None


def overlap(a, b):
    if None in (a[1], a[2], b[1], b[2]):
        return ""
    return "yes" if a[1] <= b[2] and b[1] <= a[2] else "no"


if __name__ == "__main__":
    rows = []
    for config in sorted(RUNS_DIR.glob("*/*/config.yml")):
        run_dir = config.parent
        predictor, embedder = run_dir.parent.name, run_dir.name
        out_yml = latest_out_yml(run_dir)
        if out_yml is None:
            print(f"{predictor:<26} {embedder:<12} not trained yet")
            continue
        new, old = biotrainer_metrics(out_yml), previous_metrics(predictor, embedder)
        for metric in ("accuracy", "macro_f1", "mcc"):
            n, o = new[metric], old[metric]
            rows.append([predictor, embedder, metric, f"{o[0]:.3f}", f"{n[0]:.3f}", f"{n[0] - o[0]:+.3f}",
                         f"{n[1]:.3f}-{n[2]:.3f}", f"{o[1]:.3f}-{o[2]:.3f}" if o[1] is not None else "", overlap(n, o)])

    header = ["predictor", "embedder", "metric", "previous", "biotrainer", "difference", "biotrainer_ci", "previous_ci",
              "ci_overlap"]
    OUT_TSV.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT_TSV, "w") as f:
        f.write("\t".join(header) + "\n")
        for row in rows:
            f.write("\t".join(row) + "\n")

    print(f"{'predictor':<26} {'embedder':<12} {'metric':<9} {'previous':>8} {'biotrainer':>10} {'diff':>7} {'CIs overlap':>11}")
    for row in rows:
        print(f"{row[0]:<26} {row[1]:<12} {row[2]:<9} {row[3]:>8} {row[4]:>10} {row[5]:>7} {row[8]:>11}")
    print(f"\nSaved to {OUT_TSV}")
