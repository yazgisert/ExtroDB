# train_prott5_server.py
#
# Trains a sequence-to-class classifier for each dataset split using the
# biocentral API and biotrainer (server-side training with ProtT5 embeddings).
# Unlike train_local.py, this delegates both embedding and training to the server.
#
# For each split we get test set metrics (accuracy, MCC, F1 etc.) which are
# the classification results to report in the thesis.
#

import json
import os
from biocentral_api import BiocentralAPI, SequenceTrainingData, CommonEmbedder


def parse_faa_for_training(filepath, partition):
    """
    Parse a .faa file and return a list of SequenceData objects.

    Each object carries:
      seq_id   — protein accession (e.g. WP_008296005.1)
      seq      — amino acid sequence
      label    — phenotype class (e.g. 'halophile'), parsed from the header
      set      — which partition this sequence belongs to ('train', 'val', or 'test')

    Our header format:
        >GCF_003151135.1||WP_008296005.1 description||halophile
                                                         ^^^^^^^^ this is the label
    """
    data = []
    current_id = None
    current_label = None
    current_seq = []

    with open(filepath) as f:
        for line in f:
            line = line.strip()

            if line.startswith('>'):
                # Save previous sequence
                if current_id is not None:
                    data.append(SequenceTrainingData(
                        seq_id=current_id,
                        sequence=''.join(current_seq),
                        label=current_label,
                        set=partition,
                    ))

                parts = line[1:].split('||')
                current_id = parts[1].split()[0].strip()  # protein accession
                current_label = parts[2].strip()           # phenotype label
                current_seq = []

            else:
                current_seq.append(line)

    # Don't forget the last sequence
    if current_id is not None:
        data.append(SequenceTrainingData(
            seq_id=current_id,
            sequence=''.join(current_seq),
            label=current_label,
            set=partition,
        ))

    return data


def run_training(split_name, faa_train, faa_val, faa_test, results_dir="training_results"):
    """
    Train a classifier for one split and save results to a JSON file.

    split_name  -- e.g. "split1_30pct"
    faa_*       -- paths to the three .faa files for this split
    results_dir -- folder where JSON result files are written
    """
    print(f"\nParsing sequences for {split_name} ...")
    training_data = []
    training_data += parse_faa_for_training(faa_train, "train")
    training_data += parse_faa_for_training(faa_val,   "val")
    training_data += parse_faa_for_training(faa_test,  "test")

    n_train = sum(1 for d in training_data if d.set == "train")
    n_val   = sum(1 for d in training_data if d.set == "val")
    n_test  = sum(1 for d in training_data if d.set == "test")
    labels  = sorted(set(d.label for d in training_data))
    print(f"  train: {n_train}  val: {n_val}  test: {n_test}")
    print(f"  classes ({len(labels)}): {labels}")

    print("Connecting to biocentral API ...")
    api = BiocentralAPI().wait_until_healthy()

    # sequence_to_class = multi-class classification
    # FNN = feedforward neural network (standard choice for embedding-based models)
    config = {
        "protocol":      "sequence_to_class",
        "model_choice":  "FNN",
        "embedder_name": CommonEmbedder.ProtT5,
    }

    print(f"Starting training for {split_name} (this may take a while) ...")
    result = api.train(config=config, training_data=training_data).run_with_progress()


    # Print a summary of the test results
    print(f"\n{'='*50}")
    print(f"Results for {split_name}")
    print(f"{'='*50}")
    print(result)

    # Save the full result object as JSON for later reference
    os.makedirs(results_dir, exist_ok=True)
    out_path = os.path.join(results_dir, f"{split_name}_result.json")
    try:
        with open(out_path, "w") as f:
            json.dump(result.model_dump(), f, indent=2, default=str)
        print(f"Full result saved to {out_path}")
    except Exception as e:
        print(f"Could not save result to JSON: {e}")

    return result


if __name__ == "__main__":

    # ── Size threshold test ────────────────────────────────────────────────────
    # The full 142k-sequence run fails on the server (FAILED).
    # This test uses 5,000 sequences balanced across all 5 classes (1000 per class)
    # to check if (a) multi-class training works and (b) find the size limit.
    import sys
    from collections import defaultdict

    def balanced_sample(data, n_per_class):
        by_class = defaultdict(list)
        for d in data:
            by_class[d.label].append(d)
        result = []
        for label, items in by_class.items():
            result.extend(items[:n_per_class])
        return result

    print("Loading split1_30pct for size test...")
    train_data = parse_faa_for_training("dataset_splits/rep_split1_train.faa", "train")
    val_data   = parse_faa_for_training("dataset_splits/rep_split1_val.faa",   "val")
    test_data  = parse_faa_for_training("dataset_splits/rep_split1_test.faa",  "test")

    # Take 1000 per class from train, 100 per class from val/test
    sample = balanced_sample(train_data, 1000) + \
             balanced_sample(val_data,    100) + \
             balanced_sample(test_data,   100)

    labels = sorted(set(d.label for d in sample))
    print(f"  Total: {len(sample)} sequences, classes: {labels}")

    api = BiocentralAPI().wait_until_healthy()
    config = {"protocol": "sequence_to_class", "model_choice": "FNN", "embedder_name": CommonEmbedder.ProtT5}
    result = api.train(config=config, training_data=sample).run_with_progress()
    print("Result:", result)
    sys.exit(0)
    # ── End size test ──────────────────────────────────────────────────────────

    runs = [
        ("split1_30pct",
         "dataset_splits/rep_split1_train.faa",
         "dataset_splits/rep_split1_val.faa",
         "dataset_splits/rep_split1_test.faa"),

        ("split2_30pct",
         "dataset_splits/rep_split2_train.faa",
         "dataset_splits/rep_split2_val.faa",
         "dataset_splits/rep_split2_test.faa"),

        ("split1_20pct",
         "dataset_splits/rep_split1_20pct_train.faa",
         "dataset_splits/rep_split1_20pct_val.faa",
         "dataset_splits/rep_split1_20pct_test.faa"),

        ("split2_20pct",
         "dataset_splits/rep_split2_20pct_train.faa",
         "dataset_splits/rep_split2_20pct_val.faa",
         "dataset_splits/rep_split2_20pct_test.faa"),
    ]

    for split_name, faa_train, faa_val, faa_test in runs:
        run_training(split_name, faa_train, faa_val, faa_test)
