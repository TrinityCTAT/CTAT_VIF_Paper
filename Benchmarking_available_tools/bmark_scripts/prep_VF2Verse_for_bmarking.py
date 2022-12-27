#!/usr/bin/env python

import sys, os, re
import csv


usage = f"usage: {sys.argv[0]} VF2Verse.agg.tsv\n\n"
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
    entry_counter = 0
    for row in reader:
        entry_counter += 1
        sample_name = row["sample_name"]
        virus = row["chrB"]
        virus_brkpt = row["coordB"]
        human_chr = row["chrA"]
        human_brkpt = row["coordA"]
        bmark_group = f"bmark-{virus}:{virus_brkpt}"
        read_count = row["total"]

        writer.writerow(
            {
                "sample_name": sample_name,
                "bmark_group": bmark_group,
                "human_chr": human_chr,
                "human_brkpt": human_brkpt,
                "virus": virus,
                "virus_brkpt": virus_brkpt,
                "read_count": read_count,
            }
        )

sys.exit(0)
