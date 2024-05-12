#!/usr/bin/env python3

import sys, os, re
import argparse
import pysam
from collections import defaultdict
import pyranges as pr
import csv
import subprocess
import random

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

    parser.add_argument(
        "--randomize_location",
        action="store_true",
        default=False,
        help="randomize genome location info to generate null distribution",
    )

    args = parser.parse_args()

    region_size = args.region_size
    insertions_tsv = args.insertions_tsv
    output_filename = args.output_filename
    # gistic_matrix_bdbs = args.gistic_matrix_bdbs
    cnv_regions_tsv = args.cnv_regions_tsv_tabix_gz
    randomize_location_flag = args.randomize_location

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

        if randomize_location_flag:
            (chrom, coord) = get_random_genome_pos()

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


## Random genome position code

chr_lengths = {
    "chr1": 248956422,
    "chr2": 242193529,
    "chr3": 198295559,
    "chr4": 190214555,
    "chr5": 181538259,
    "chr6": 170805979,
    "chr7": 159345973,
    "chr8": 145138636,
    "chr9": 138394717,
    "chr10": 133797422,
    "chr11": 135086622,
    "chr12": 133275309,
    "chr13": 114364328,
    "chr14": 107043718,
    "chr15": 101991189,
    "chr16": 90338345,
    "chr17": 83257441,
    "chr18": 80373285,
    "chr19": 58617616,
    "chr20": 64444167,
    "chr21": 46709983,
    "chr22": 50818468,
    "chrX": 156040895,
    "chrY": 57227415,
}

chromosomes = list(chr_lengths.keys())

sum_genome_length = sum(chr_lengths.values())

ordered_chromosomes = list()
cumsum = 0
for chromosome in chromosomes:
    chr_len = chr_lengths[chromosome]
    ordered_chromosomes.append(
        {"chrom": chromosome, "lend": cumsum + 1, "rend": cumsum + chr_len}
    )
    cumsum += chr_len


def get_random_genome_pos():

    rand_chrom_pos = random.randint(1, sum_genome_length)

    for chrom_struct in ordered_chromosomes:
        if (
            rand_chrom_pos >= chrom_struct["lend"]
            and rand_chrom_pos <= chrom_struct["rend"]
        ):
            chrom = chrom_struct["chrom"]
            chrom_pos = rand_chrom_pos - chrom_struct["lend"] + 1

            return (chrom, chrom_pos)

    raise RuntimeError("Error, no random position selected - shouldn't happen - bug!")


if __name__ == "__main__":
    main()
