# Copies the ProtT5 algae results from eukaryote_embedder_clustered_result.json into the file format the algae plot
# scripts read (filtered and unfiltered).

import json

SOURCE = "results/training_results/eukaryote_embedder_clustered_result.json"
OUT    = {"filtered":   "results/training_results/binary_eukaryote_eval_clustered_result.json",
          "unfiltered": "results/training_results/binary_eukaryote_eval_unfiltered_clustered_result.json"}


def legacy_format(result, n):
    # Same structure as binary_eukaryote_eval_filtered_result.json, so the existing plot scripts can read it
    def entry(r):
        return {"n_proteins": r["n_proteins"], "recall": r["recall"], "mean_prob": r["mean_prob"],
                "recall_ci_lower": r["ci_lower"], "recall_ci_upper": r["ci_upper"]}
    return {
        "task": "binary ProtT5 model (clustered split) on eukaryote algae",
        "model": "results/training_results/binary_prott5_clustered_model.pt",
        "n_algae_proteins": n,
        "overall": entry(result["overall"]),
        "by_class": {c: entry(r) for c, r in result["by_class"].items()},
        "by_species": result["by_species"],
    }


if __name__ == "__main__":
    with open(SOURCE) as f:
        data = json.load(f)
    prott5 = data["results"]["prott5"]
    for key, path in OUT.items():
        result = prott5 if key == "filtered" else prott5["unfiltered"]
        n = data["n_algae_proteins"] if key == "filtered" else prott5["unfiltered"]["n_algae_proteins"]
        with open(path, "w") as f:
            json.dump(legacy_format(result, n), f, indent=2)
        print(f"Saved to {path}")
