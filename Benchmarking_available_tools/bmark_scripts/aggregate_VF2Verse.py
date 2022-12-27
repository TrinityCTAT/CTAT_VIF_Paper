#!/usr/bin/env python3

import pandas as pd
import pysam
import os, sys, re
import logging
import sys, time
import itertools

logging.basicConfig(format="\n %(levelname)s : %(message)s", level=logging.DEBUG)
logger = logging.getLogger(__name__)


def reformat(input_file, sample_name, virus_type):

    ## read in the input file
    df = pd.read_csv(input_file, sep="\t")

    header = [
        "chrA",
        "coordA",
        "orientA",
        "chrB",
        "coordB",
        "orientB",
        "Support_reads_pair_softclip",
        "Confidence",
    ]
    df.columns = header
    df

    ## Create sample name
    df.insert(0, "sample_name", sample_name)

    ## Create contig Name
    ### adjust orientB if not given
    idx = df["orientB"].isnull()
    df.loc[idx, ["orientB"]] = "NA"
    ### combine
    contig = (
        df["chrA"]
        + "~"
        + df["coordA"].astype(str)
        + "~"
        + df["orientA"]
        + "~"
        + df["chrB"]
        + "~"
        + df["coordB"].astype(str)
        + "~"
        + df["orientB"]
    )
    df.insert(1, "entry", contig)

    ## get the total number of suporting  reads
    sums = (
        df["Support_reads_pair_softclip"]
        .str.split("+", expand=True)
        .astype(int)
        .sum(axis=1)
    )
    df["total"] = sums

    idx = df.chrA.str.match("^chrVirus")
    df.loc[idx, ["chrA", "coordA", "orientA", "chrB", "coordB", "orientB"]] = df.loc[
        idx, ["chrB", "coordB", "orientB", "chrA", "coordA", "orientA"]
    ].values

    df["chrB"] = df["chrB"].replace({"chrVirus": virus_type})

    return df


def main():

    usage = f"usage: {sys.argv[0]} fileA [fileB ...]\n\n"
    if len(sys.argv) < 2:
        exit(usage)

    files = sys.argv[1:]

    concat_df = None
    for file in files:
        sample_name = os.path.basename(file).replace("_results-virus-loci.txt", "")
        df = reformat(input_file=file, sample_name=sample_name, virus_type=sample_name)
        concat_df = pd.concat([concat_df, df]) if concat_df is not None else df

    concat_df.to_csv(sys.stdout, sep="\t", index=False)

    sys.exit(0)


if __name__ == "__main__":
    main()
