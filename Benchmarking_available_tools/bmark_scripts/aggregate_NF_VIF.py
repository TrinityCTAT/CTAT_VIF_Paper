#!/usr/bin/env python3

import pandas as pd
import pysam
import os, sys, re
import logging
import sys, time
import itertools
import argparse


logging.basicConfig(format="\n %(levelname)s : %(message)s", level=logging.DEBUG)
logger = logging.getLogger(__name__)


def readFiles(dirname):

    # Get the files in the output directory
    file_list = os.listdir(dirname)
    # select the filtered outputs for each virus
    # r = re.compile(".*_table_filtered.csv")
    r = re.compile(".*_table.csv")
    newlist = list(filter(r.match, file_list))
    df = pd.DataFrame()
    for table_file in newlist:
        output_file = os.path.join(dirname, table_file)
        tmp = pd.read_csv(output_file)
        df = pd.concat([df, tmp])

    return df


def main():

    usage = f"\n\n\tusage: {sys.argv[0]} nf_vif_outdir [nf_vif_outdir, ...]\n\n"
    if len(sys.argv) < 2:
        exit(usage)

    concat_df = pd.DataFrame()
    outdirs = sys.argv[1:]
    for outdir in outdirs:
        df = readFiles(outdir)
        sample_name = os.path.basename(outdir).replace(".nf-vif", "")
        df["sample_name"] = sample_name
        concat_df = pd.concat([concat_df, df])

    concat_df.to_csv(sys.stdout, sep="\t", index=False)

    sys.exit(0)


if __name__ == "__main__":
    main()
