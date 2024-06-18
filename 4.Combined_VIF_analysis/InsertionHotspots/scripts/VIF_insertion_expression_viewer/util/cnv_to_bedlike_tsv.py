#!/usr/bin/env python

import pandas as pd
import sys, os, re
import argparse


parser = argparse.ArgumentParser(
    description="convert set of cnv tsvs into tabix bed",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
)

parser.add_argument(
    "--cnv_tsv", type=str, required=True, nargs="+", help="cnv tsv file"
)


parser.add_argument(
    "--output", type=str, required=True, help="output filename (ie. cnvs.bed)"
)


args = parser.parse_args()


tsv_files = args.cnv_tsv
output_filename = args.output

data = None

for tsv_file in tsv_files:

    tsv_data = pd.read_csv(tsv_file, sep="\t")
    if data is None:
        data = tsv_data
    else:
        data = pd.concat([data, tsv_data])

# reorder columns
data = data[["Chrom", "Start", "End", "sample", "value"]]

data.sort_values(["Chrom", "Start", "End"], inplace=True)

data.to_csv(output_filename, sep="\t", header=False, index=False)

sys.exit(0)
