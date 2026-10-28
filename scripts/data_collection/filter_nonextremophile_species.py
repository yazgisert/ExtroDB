# Removes the non-extremophile proteins of the species flagged by flag_nonextremophile_species.py. The result is the
# input for cluster_nonextremophiles.py.

import csv
import re

from flag_nonextremophile_species import OUT_TSV as FLAGS_TSV, species_name

NONEXTREM_FAA = "nonextremophile_data/nonextremophile_filtered.faa"
OUT_FAA       = "nonextremophile_data/nonextremophile_species_filtered.faa"


if __name__ == "__main__":
    with open(FLAGS_TSV) as f:
        flagged = {row["species"] for row in csv.DictReader(f, delimiter="\t") if row["flags"]}

    n_kept = n_removed = 0
    with open(NONEXTREM_FAA) as f, open(OUT_FAA, "w") as out:
        keep = True
        for line in f:
            if line.startswith(">"):
                match = re.search(r" OS=(.+?) OX=", line)
                keep = not (match and species_name(match.group(1)) in flagged)
                n_kept += keep
                n_removed += not keep
            if keep:
                out.write(line)

    print(f"{len(flagged)} flagged species: removed {n_removed:,} proteins, kept {n_kept:,}")
    print(f"Saved to {OUT_FAA}")
