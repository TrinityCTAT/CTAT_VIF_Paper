#!/usr/bin/env python3

import sys, os, re
import pandas as pd


## See: https://github.com/namphuon/ViFi
## for formatting documentation


def ReformatVifi(input_file, sample_tag):

    with open(input_file, "rt") as fh:

        df = pd.read_csv(
            fh,
            comment="#",
            sep="\t",
            header=None,
            names=[
                "chr",
                "minpos",
                "maxpos",
                "reads",
                "forward",
                "reverse",
                "Min",
                "Max",
                "Split1",
                "Split2",
            ],
        )
        df.insert(0, "sample_name", sample_tag)

    return df


def main():

    usage = f"usage: {sys.argv[0]} fileA [fileB ...]\n\n"
    if len(sys.argv) < 2:
        exit(usage)

    files = sys.argv[1:]

    # initiate the ViFi object
    concat_df = None
    for file in files:
        sample_name = os.path.basename(file).replace(".clusters.txt", "")
        df = ReformatVifi(input_file=file, sample_tag=sample_name)
        concat_df = pd.concat([concat_df, df]) if concat_df is not None else df

    concat_df.to_csv(sys.stdout, sep="\t", index=False)

    sys.exit(0)


if __name__ == "__main__":
    main()
