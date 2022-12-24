#!/usr/bin/env python3

import os, sys, re
import pandas as pd


def reformat(input_file, sample_tag):

    message_str = f"\t\tReformating the file."
    print(message_str, file=sys.stderr)

    #####################
    # Put into dataframe
    #####################
    # column_names = ["chrA","Min","Max","total_reads", "L_reads","R_reads", "min","max","split1","split2"]
    df = pd.read_csv(input_file, sep="\t", index_col=False)

    df.insert(1, "#sample", sample_tag)
    df = df.drop(labels=["LIB"], axis=1)

    # Rename the columns
    df = df.rename(
        {
            "Chr": "chrA",
            "Human Pos": "coordA",
            "Sign": "orientA",
            "Viral Sign": "orientB",
            "Viral Pos": "coordB",
            "Read Count": "total",
            # replace spaceswith underscores
            "Split Reads": "Split_Reads",
            "Uniquely Mapped Reads": "Uniquely_Mapped_Reads",
            "Multiply Mapped Reads": "Multiply_Mapped_Reads",
            "Rank1 Hits": "Rank1_Hits",
        },
        axis="columns",
    )

    # chrA : adjust the format, add 'chr' to begining
    df["chrA"] = "chr" + df["chrA"].astype(str)

    # Add virus of interest

    df.insert(5, "chrB", sample_tag)

    # entry : create and add the entry column
    entry = (
        df["chrA"]
        + "~"
        + df["coordA"].astype("str")
        + "~"
        + df["chrB"]
        + "~"
        + df["coordB"].astype("str")
    )
    df.insert(1, "entry", entry)

    # Reorder columns
    column_order = [
        "#sample",
        "entry",
        "chrA",
        "coordA",
        "orientA",
        "chrB",
        "coordB",
        "orientB",
        "total",
        "Split_Reads",
        "Uniquely_Mapped_Reads",
        "Multiply_Mapped_Reads",
        "Rank1_Hits",
    ]
    df = df[column_order]

    return df


def main():

    usage = f"usage: {sys.argv[0]} fileA [fileB ...]\n\n"
    if len(sys.argv) < 2:
        exit(usage)

    concat_df = None
    files = sys.argv[1:]
    for filename in files:
        sample_name = os.path.basename(filename).replace(".final_hits.txt", "")
        df = reformat(input_file=filename, sample_tag=sample_name)
        concat_df = pd.concat([concat_df, df]) if concat_df is not None else df

    concat_df.to_csv(sys.stdout, sep="\t", index=False)

    sys.exit(0)


if __name__ == "__main__":
    main()
