#!/usr/bin/env python

import sys, os, re
import csv


usage = f"usage: {sys.argv[0]} vifi.agg.tsv\n\n"
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
        virus = sample_name
        bmark_group = f"bmark-{virus}-entry{entry_counter}"
        read_count = row["reads"]

        brkpts = []
        if int(row["Split1"]) > 0:
            brkpts.append(row["Split1"])
        if int(row["Split2"]) > 0:
            brkpts.append(row["Split2"])
        if len(brkpts) < 1:
            brkpts.append(int((int(row["minpos"]) + int(row["maxpos"])) / 2))

        for brkpt in brkpts:
            writer.writerow(
                {
                    "sample_name": sample_name,
                    "bmark_group": bmark_group,
                    "human_chr": row["chr"],
                    "human_brkpt": brkpt,
                    "virus": virus,
                    "virus_brkpt": 1,
                    "read_count": read_count,
                }
            )

sys.exit(0)
