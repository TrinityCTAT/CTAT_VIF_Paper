#!/usr/bin/env python3

import sys, os, re
import argparse
import pysam
from collections import defaultdict
import csv
import subprocess
import glob
import statistics
from statsmodels.distributions.empirical_distribution import ECDF

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
        "--expr_bed_tabix_dir",
        type=str,
        required=True,
        help="dir containing the gene expr beds bgzip'd and tabix indexed",
    )

    args = parser.parse_args()

    region_size = args.region_size
    insertions_tsv = args.insertions_tsv
    output_filename = args.output_filename
    expr_bed_tabix_dir = args.expr_bed_tabix_dir

    tcga_to_expr_bed_filename = get_TCGA_bed_file_mappings(expr_bed_tabix_dir)

    insertions_tsv_fh = open(insertions_tsv, "rt")
    tab_reader = csv.DictReader(insertions_tsv_fh, delimiter="\t")

    expr_info_ofh = open(output_filename, "wt")
    tab_writer = csv.DictWriter(
        expr_info_ofh,
        fieldnames=[
            "TCGA",
            "sample_id",
            "chrom",
            "virus",
            "contig",
            "sum_expr",
            "quantile",
            "fold_change",
        ],
        delimiter="\t",
    )
    tab_writer.writeheader()

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

        if tcga not in tcga_to_expr_bed_filename:
            continue

        expr_bed_filename = tcga_to_expr_bed_filename[tcga]

        sample_names_filename = expr_bed_filename.replace(
            "gene_expr_gene_span_list.sorted.gz", "sample_names"
        )
        sample_list = None
        with open(sample_names_filename, "rt") as fh:
            line = next(fh)
            line = line.rstrip()
            sample_list = line.split(",")

        print(sample_list)
        for i, sample in enumerate(sample_list):
            print("\t".join([str(i), sample]))

        expr_regions_tabix = pysam.TabixFile(expr_bed_filename)

        expr_sums = defaultdict(int)
        for expr_entry in expr_regions_tabix.fetch(chrom, region_lend, region_rend):
            print(expr_entry)
            (
                expr_region_chrom,
                expr_region_lend,
                expr_region_rend,
                gene,
                sample_expr_vals,
            ) = expr_entry.split("\t")

            sample_expr_vals = sample_expr_vals.split(",")
            print(sample_expr_vals)

            assert len(sample_expr_vals) == len(
                sample_list
            ), "Error, {} size of sample expr vals != num samples {}".format(
                len(sample_expr_vals), len(sample_list)
            )

            sample_expr_vals = [float(i) for i in sample_expr_vals]

            for i, sample in enumerate(sample_list):
                expr_sums[sample] += sample_expr_vals[i]

        sample_expr_val = expr_sums[sample]

        all_expr_vals = list(expr_sums.values())
        mean_expr_val = statistics.mean(all_expr_vals)

        ecdf = ECDF(all_expr_vals)
        # sample_expr_val = max(sample_expr_val, 0.001)

        sample_q = ecdf(sample_expr_val)

        print("sample_expr_val: {}".format(sample_expr_val))
        print("mean_expr_val: {}".format(mean_expr_val))
        print("sample_q: {}".format(sample_q))

        """

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
        """

    expr_info_ofh.close()
    sys.exit(0)


def get_TCGA_bed_file_mappings(expr_bed_dir):

    files = glob.glob(expr_bed_dir + "/TCGA*.genesym*.tsv.gz*tbi")

    tcga_bed_file_mappings = dict()
    for filename in files:
        filename = filename.replace(".tbi", "")
        bname = os.path.basename(filename)
        tcga_tok = bname.split(".")[0].split("-")[1]
        tcga_bed_file_mappings[tcga_tok] = filename

    return tcga_bed_file_mappings


if __name__ == "__main__":
    main()
