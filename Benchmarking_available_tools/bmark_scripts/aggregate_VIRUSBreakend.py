#!/usr/bin/env python3

import pandas as pd
import pysam
import os, sys, re
import logging
import sys, time

logging.basicConfig(format="\n %(levelname)s : %(message)s", level=logging.DEBUG)
logger = logging.getLogger(__name__)


def readFile(input_file):
    file_open = open(input_file, "r")
    file_input = file_open.readlines()
    file_open.close()

    new = []
    for i in file_input:
        if i[0:2] != "##":
            new.append(i.rstrip().split("\t"))
    # get the header
    header = new.pop(0)
    # create new dataframe
    df = pd.DataFrame(new)
    df.columns = header

    return df


def reformatOutput(input_file, virus_type):

    integrations = readFile(input_file)

    # Fix the ID column for merging
    integrations["ID"] = integrations["ID"].str.replace("_host", "")
    integrations["ID"] = integrations["ID"].str.replace("_virus", "")

    # seperate the human and virus integrations
    # chr_idx = self.integrations["#CHROM"].str.contains("chr")
    chr_idx = integrations["REF"].str.contains("N")
    human_idx = [i for i, j in enumerate(chr_idx) if j == True]
    viral_idx = [i for i, j in enumerate(chr_idx) if j == False]

    human_df = integrations.iloc[human_idx]

    # combine the viral POS to the correct human integration point
    viral_df = integrations.iloc[viral_idx]
    viral_df = viral_df.rename({"POS": "coordB"}, axis="columns")
    viral_df_pos = viral_df[["coordB", "ID"]]
    df = human_df.merge(viral_df_pos, on="ID")

    coordB = df.pop("coordB")
    df.insert(3, "coordB", coordB)

    # Adjust column names
    df = df.rename({"#CHROM": "chrA", "POS": "coordA"}, axis="columns")

    # human orientation
    orientation = [i.split("|")[1] for i in df["INFO"]]
    df.insert(3, "orientA", orientation)
    # insert virus type
    df.insert(4, "chrB", virus_type)

    # ~~~~~~~~~~~~~~
    # Get total reads
    # ~~~~~~~~~~~~~~
    BVF_list = []
    r = re.compile("BVF=*")
    for i in viral_df["INFO"]:
        tmp = i.split(";")
        BVF = list(filter(r.match, tmp))[0]
        BVF = BVF.replace("BVF=", "")
        BVF_list.append(BVF)

    df.insert(6, "total", BVF_list)

    # ~~~~~~~~~~
    # Create entry
    # ~~~~~~~~~~
    entry = (
        df["chrA"]
        + "~"
        + df["coordA"].astype(str)
        + "~"
        + df["orientA"]
        + "~"
        + df["chrB"]
        + "~"
        + df["coordB"].astype(str)
    )
    df.insert(0, "entry", entry)

    # ~~~~~~~~~~
    ## Create sample name
    # ~~~~~~~~~~
    # df.insert(0, "#sample", [self.sample_tag] * df.shape[0])

    df.insert(0, "sample_name", virus_type)

    df = df.drop(labels=["ID"], axis=1)

    return df


def main():

    usage = f"usage: {sys.argv[0]} fileA [fileB ...]\n\n"
    if len(sys.argv) < 2:
        exit(usage)

    concat_df = None
    files = sys.argv[1:]
    for filename in files:
        sample_name = os.path.basename(filename).replace(".virusbreakend.vcf", "")
        df = reformatOutput(filename, sample_name)
        concat_df = pd.concat([concat_df, df]) if concat_df is not None else df

    concat_df.to_csv(sys.stdout, sep="\t", index=False)

    sys.exit(0)


if __name__ == "__main__":
    main()
