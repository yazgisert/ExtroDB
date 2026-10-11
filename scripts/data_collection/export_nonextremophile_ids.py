# Writes the exact non-extremophile UniProt sample (IDs + headers) with flags for the extremophile filter, the
# embeddings and the 30% cluster representatives. The sample can't be redrawn later because UniProt changes.

import h5py

RAW_FAA      = "nonextremophile_data/nonextremophile.faa"
FILTERED_FAA = "nonextremophile_data/nonextremophile_filtered.faa"
EMBEDDINGS   = "embeddings/nonextremophile.h5"
CLUSTERS     = "nonextremophile_data/clustering_30pct/nonext_cluster_cluster.tsv"
OUT          = "curation/source_tables/nonextremophile_uniprot_sample.tsv"


def headers(path):
    with open(path) as f:
        return [line[1:].strip() for line in f if line.startswith(">")]


if __name__ == "__main__":
    filtered = {h.split()[0] for h in headers(FILTERED_FAA)}
    with h5py.File(EMBEDDINGS, "r") as f:
        embedded = set(f.keys())
    with open(CLUSTERS) as f:
        representatives = {line.split("\t")[0] for line in f}

    with open(OUT, "w") as out:
        out.write("protein_id\tuniprot_header\tkept_after_extremophile_filter\tembedded\tcluster_representative_30pct\n")
        for header in headers(RAW_FAA):
            protein_id = header.split()[0]
            description = header.split(" ", 1)[1].split("||")[0] if " " in header else ""
            out.write(f"{protein_id}\t{description}\t{protein_id in filtered}\t{protein_id in embedded}\t"
                      f"{protein_id in representatives}\n")
    print(f"Saved to {OUT}")
