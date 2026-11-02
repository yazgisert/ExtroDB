# Removes the non-extremophile proteins of the species flagged by flag_nonextremophile_species.py, and of the species
# without any growth data whose genus contains one of our prokaryotic extremophiles. The result is the input for
# cluster_nonextremophiles.py.

import csv
import re

from flag_nonextremophile_species import OUT_TSV as FLAGS_TSV, species_name

NONEXTREM_FAA = "nonextremophile_data/nonextremophile_filtered.faa"
MASTER_INDEX  = "searchable_resource/extremophile_master_index.tsv"
OUT_FAA       = "nonextremophile_data/nonextremophile_species_filtered.faa"


def extremophile_genera():
    with open(MASTER_INDEX) as f:
        return {row["taxon_name"].replace('"', "").split()[0]
                for row in csv.DictReader(f, delimiter="\t") if row["domain_scope"] == "prokaryote"}


def has_growth_data(row):
    return any(row[k] for k in ("temperature_optima", "ph_optima", "nacl_optima_M", "tempura_ogt"))


if __name__ == "__main__":
    genera = extremophile_genera()
    with open(FLAGS_TSV) as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    flagged = {r["species"] for r in rows if r["flags"]}
    no_data = {r["species"] for r in rows
               if not r["flags"] and r["species"].split()[0] in genera and not has_growth_data(r)}
    removed = flagged | no_data

    n_kept = n_removed = 0
    with open(NONEXTREM_FAA) as f, open(OUT_FAA, "w") as out:
        keep = True
        for line in f:
            if line.startswith(">"):
                match = re.search(r" OS=(.+?) OX=", line)
                keep = not (match and species_name(match.group(1)) in removed)
                n_kept += keep
                n_removed += not keep
            if keep:
                out.write(line)

    print(f"{len(flagged)} flagged species and {len(no_data)} species of extremophile genera without growth data: "
          f"removed {n_removed:,} proteins, kept {n_kept:,}")
    print(f"Saved to {OUT_FAA}")
