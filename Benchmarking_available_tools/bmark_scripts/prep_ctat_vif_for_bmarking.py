#!/usr/bin/env python

import sys, os, re
import csv


usage = f"usage: {sys.argv[0]} vif.agg.tsv\n\n"
if len(sys.argv) < 2:
    exit(usage)


out_colnames = [
    "sample_name",
    "human_chr",
    "human_brkpt",
    "virus",
    "virus_brkpt",
    "bmark_group",
    "read_count",
]

writer = csv.DictWriter(sys.stdout, fieldnames=out_colnames, delimiter="\t")
writer.writeheader()

fname = sys.argv[1]
with open(fname, "rt") as fh:

    reader = csv.DictReader(fh, delimiter="\t")

    for row in reader:
        sample_name = row["sample_name"]
        bmark_group = row["virus_brkend_grp"]
        virus = bmark_group.split(":")[0]
        (human_chr, human_coord, virus_coord) = (
            (row["chrA"], row["coordA"], row["coordB"])
            if row["chrB"] == virus
            else (row["chrB"], row["coordB"], row["coordA"])
        )
        read_count = row["total"]

        writer.writerow(
            {
                "sample_name": sample_name,
                "bmark_group": bmark_group,
                "human_chr": human_chr,
                "human_brkpt": human_coord,
                "virus": virus,
                "virus_brkpt": virus_coord,
                "read_count": read_count,
            }
        )

sys.exit(0)
