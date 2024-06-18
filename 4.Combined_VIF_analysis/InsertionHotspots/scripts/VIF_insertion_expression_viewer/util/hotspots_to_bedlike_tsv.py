#!/usr/bin/env python

import pandas as pd
import sys, os, re


usage = f"usage: {sys.argv[0]} hotspot_sample_data.refined.sample_data.tsv"

if len(sys.argv) < 2:
    exit(usage)

tsv_file = sys.argv[1]

data = pd.read_csv(tsv_file, sep="\t")


def extract_insertion_info(val):
    # chr8~127979782~+~HPV18~462~-

    chrA, coordA, orientA, chrB, coordB, orientB = val.split("~")

    if re.match("^chr", chrA):

        (human_chr, human_coord, human_orient, virus, virus_coord, virus_orient) = (
            chrA,
            coordA,
            orientA,
            chrB,
            coordB,
            orientB,
        )

    else:
        (virus, virus_coord, virus_orient, human_chr, human_coord, human_orient) = (
            chrA,
            coordA,
            orientA,
            chrB,
            coordB,
            orientB,
        )

    return pd.Series(
        {
            "chr": human_chr,
            "start": int(human_coord),
            "end": int(human_coord) + 1,
            "human_orient": human_orient,
            "virus": virus,
            "virus_coord": virus_coord,
            "virus_orient": virus_orient,
        }
    )


extracted_coords_n_orients_df = data["contig"].apply(extract_insertion_info)

merged_df = pd.concat(
    [extracted_coords_n_orients_df, data[["total", "sample_id", "contig"]]], axis=1
)

outfilename = tsv_file + ".bedlike.tsv"

merged_df.sort_values(["chr", "start", "end"], inplace=True)

merged_df.to_csv(outfilename, sep="\t", header=False, index=False)
