#!/usr/bin/env python3

import sys, os, re
import pandas as pd


## See: https://github.com/namphuon/ViFi
## for formatting documentation


def ReformatFastVifi(input_file, range_file, sample_tag):

    with open(input_file, "rt") as fh:

        df = pd.read_csv(
            fh,
            comment="#",
            sep="\t",
            header=None,
            names=["chr", "minpos", "maxpos", "reads", "forward", "reverse"],
        )
        df.insert(0, "sample_name", sample_tag)

    with open(range_file, "rt") as fh:
        range_df = pd.read_csv(fh, sep=",")
        range_df = range_df.rename({"Chr": "chrom"})
        # range_df = range_df.drop([0], axis=0)

    df = pd.concat([df.reset_index(drop=True), range_df.reset_index(drop=True)], axis=1)

    return df


def main():

    usage = f"usage: {sys.argv[0]} fileA.clusters.txt [fileB.clusters.txt ...]\n\n"
    if len(sys.argv) < 2:
        exit(usage)

    files = sys.argv[1:]

    # initiate the ViFi object
    concat_df = None
    for file in files:
        range_file = file + ".range"
        sample_name = os.path.basename(os.path.dirname(file)).replace("fastvifi.", "")
        df = ReformatFastVifi(
            input_file=file, range_file=range_file, sample_tag=sample_name
        )
        concat_df = pd.concat([concat_df, df]) if concat_df is not None else df

    concat_df.to_csv(sys.stdout, sep="\t", index=False)

    sys.exit(0)


if __name__ == "__main__":
    main()
