#!/usr/bin/env python3

import sys, os, re
import argparse
import pysam
from collections import defaultdict
import pyranges as pr
import csv
import subprocess

UTILDIR = os.path.join(os.path.dirname(__file__), "util")


def main():

    parser = argparse.ArgumentParser(
        description="extract most extreme CNV within region",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "--region_size",
        type=int,
        required=True,
        help="size of region around insertion site",
    )
    parser.add_argument(
        "--insertions_tsv", type=str, required=True, help="virus insertions tsv",
    )
    parser.add_argument(
        "--output_filename",
        type=str,
        required=True,
        help="name of output file (can include pdf or svg extension",
    )

    parser.add_argument(
        "--cnv_regions_tsv_tabix_gz",
        type=str,
        required=True,
        help="cnv regions in tsv bgzip'd and tabix indexed",
    )

    args = parser.parse_args()

    region_size = args.region_size
    insertions_tsv = args.insertions_tsv
    output_filename = args.output_filename
    # gistic_matrix_bdbs = args.gistic_matrix_bdbs
    cnv_regions_tsv = args.cnv_regions_tsv_tabix_gz

    # cnv info here: https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/CNV_Pipeline/
    # value: The GDC further transforms these copy number values into segment mean values, which are equal to log2(copy-number/ 2). Diploid regions will have a segment mean of zero, amplified regions will have positive values, and deletions will have negative values.

    insertions_tsv_fh = open(insertions_tsv, "rt")
    tab_reader = csv.DictReader(insertions_tsv_fh, delimiter="\t")

    cnv_tsv_ofh = open(output_filename, "wt")
    tab_writer = csv.DictWriter(
        cnv_tsv_ofh,
        fieldnames=[
            "TCGA",
            "sample",
            "chrom",
            "virus",
            "insertion",
            "cnv_lend",
            "cnv_rend",
            "copy_number",
        ],
        delimiter="\t",
    )
    tab_writer.writeheader()

    cnv_regions_tabix = pysam.TabixFile(cnv_regions_tsv)

    for row in tab_reader:

        sample = row["sample_id"]
        cohort = row["cohort"]
        if cohort != "TCGA":
            continue

        tcga = row["project"]
        chrom = row["humanchr"]
        coord = int(row["human_coord"])
        virus = row["virus"]
        insertion = row["contig"]

        region_lend, region_rend = max(1, coord - region_size), coord + region_size

        largest_cnv = None

        for cnv_region in cnv_regions_tabix.fetch(chrom, region_lend, region_rend):

            (
                cnv_region_chrom,
                cnv_region_lend,
                cnv_region_rend,
                sample_name,
                value,
            ) = cnv_region.split("\t")

            if not sample_name == sample:
                continue

            cnv_region_lend = int(cnv_region_lend)
            cnv_region_rend = int(cnv_region_rend)

            value = float(value)

            copy_number = 2 * 2 ** value

            if largest_cnv is None or largest_cnv["copy_number"] < copy_number:

                cnv_region_dict = {
                    "chrom": chrom,
                    "region_lend": cnv_region_lend,
                    "region_rend": cnv_region_rend,
                    "copy_number": copy_number,
                }

                largest_cnv = cnv_region_dict

        if largest_cnv is not None:

            tab_writer.writerow(
                {
                    "TCGA": tcga,
                    "sample": sample,
                    "chrom": chrom,
                    "virus": virus,
                    "insertion": insertion,
                    "cnv_lend": largest_cnv["region_lend"],
                    "cnv_rend": largest_cnv["region_rend"],
                    "copy_number": largest_cnv["copy_number"],
                }
            )

    cnv_tsv_ofh.close()
    sys.exit(0)


if __name__ == "__main__":
    main()
