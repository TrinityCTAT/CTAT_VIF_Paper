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
import numpy
import math
import random
import logging


logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s : %(levelname)s : %(message)s",
    datefmt="%H:%M:%S",
)
logger = logging.getLogger(__name__)


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

    parser.add_argument(
        "--by_gene_in_region",
        action="store_true",
        default=False,
        help="perform per-gene analysis within region instead of cumulatively (default)",
    )

    parser.add_argument(
        "--randomize_location",
        action="store_true",
        default=False,
        help="randomize genome location info to generate null distribution",
    )

    parser.add_argument(
        "--randomize_sample",
        action="store_true",
        default=False,
        help="randomize target sample_id to generate null distribution",
    )

    parser.add_argument(
        "--rand_iters",
        type=int,
        default=10,
        help="number of random iterations per entry",
    )

    ## parse args
    args = parser.parse_args()

    region_size = args.region_size
    insertions_tsv = args.insertions_tsv
    output_filename = args.output_filename
    expr_bed_tabix_dir = args.expr_bed_tabix_dir
    by_gene_flag = args.by_gene_in_region
    randomize_location_flag = args.randomize_location
    randomize_sample_flag = args.randomize_sample
    rand_iters = args.rand_iters

    ## do work:
    project_to_expr_bed_filename = get_project_bed_file_mappings(expr_bed_tabix_dir)

    insertions_tsv_fh = open(insertions_tsv, "rt")
    tab_reader = csv.DictReader(insertions_tsv_fh, delimiter="\t")

    expr_info_ofh = open(output_filename, "wt")
    tab_writer = csv.DictWriter(
        expr_info_ofh,
        fieldnames=[
            "cohort",
            "project",
            "sample_id",
            "seqtype",
            "chrom",
            "virus",
            "contig",
            "sample_region_expr",
            "mean_region_expr",
            "log2_fold_change",
            "expr_quantile",
            "region_gene_list",
        ],
        delimiter="\t",
    )
    tab_writer.writeheader()

    for row in tab_reader:

        sample = row["sample_id"]
        cohort = row["cohort"]

        if cohort not in ["TCGA", "GTEx"]:
            # only tcga and gtex have expression matrices.
            continue

        project = row["project"]

        if project not in project_to_expr_bed_filename:
            logger.error(f"-missing {project} expr bed filename")
            continue

        iters = 1
        if randomize_sample_flag or randomize_location_flag:
            iters = rand_iters

        for i in range(iters):
            examine_row(row, args, tab_writer, project_to_expr_bed_filename)

    expr_info_ofh.close()

    sys.exit(0)


def examine_row(row, args, tab_writer, project_to_expr_bed_filename):

    region_size = args.region_size
    expr_bed_tabix_dir = args.expr_bed_tabix_dir
    by_gene_flag = args.by_gene_in_region
    randomize_location_flag = args.randomize_location
    randomize_sample_flag = args.randomize_sample

    sample = row["sample_id"]
    chrom = row["humanchr"]
    coord = int(row["human_coord"])
    virus = row["virus"]
    insertion = row["contig"]
    project = row["project"]
    cohort = row["cohort"]

    if randomize_location_flag:
        (chrom, coord) = get_random_genome_pos()

    region_lend, region_rend = max(1, coord - region_size), coord + region_size

    expr_bed_filename = project_to_expr_bed_filename[project]
    # print("-searching {}".format(expr_bed_filename))

    sample_names_filename = expr_bed_filename.replace(
        "gene_expr_gene_span_list.sorted.gz", "sample_names"
    )
    sample_list = None
    with open(sample_names_filename, "rt") as fh:
        line = next(fh)
        line = line.rstrip()
        sample_list = line.split(",")

    if randomize_sample_flag:
        sample = random.choice(sample_list)

    # print(sample_list)
    # for i, sample in enumerate(sample_list):
    #    print("\t".join([str(i), sample]))

    expr_regions_tabix = pysam.TabixFile(expr_bed_filename)

    expr_sums = defaultdict(int)
    gene_list = list()
    for expr_entry in expr_regions_tabix.fetch(chrom, region_lend, region_rend):
        # print(expr_entry)
        (
            expr_region_chrom,
            expr_region_lend,
            expr_region_rend,
            gene,
            sample_expr_vals,
        ) = expr_entry.split("\t")
        gene_list.append(gene)
        #
        if by_gene_flag:
            gene_list = [gene]
            expr_sums.clear()

        sample_expr_vals = sample_expr_vals.split(",")
        # print(sample_expr_vals)

        assert len(sample_expr_vals) == len(
            sample_list
        ), "Error, {} size of sample expr vals != num samples {}".format(
            len(sample_expr_vals), len(sample_list)
        )

        sample_expr_vals = [float(i) for i in sample_expr_vals]

        for i, sample_entry in enumerate(sample_list):
            expr_sums[sample_entry] += sample_expr_vals[i]

        if by_gene_flag:
            run_stats(expr_sums, sample, tab_writer, gene_list, row)

    if not by_gene_flag:
        run_stats(expr_sums, sample, tab_writer, gene_list, row)

    return


def run_stats(expr_sums, sample, tab_writer, gene_list, row):

    # print("len expr sums = {}".format(len(expr_sums)))
    if len(expr_sums) == 0:
        # nothing to report.
        return

    ## Examine sample expression stats
    sample_expr_val = expr_sums[sample]

    all_expr_vals = list(expr_sums.values())
    mean_expr_val = statistics.mean(all_expr_vals)

    all_expr_vals_minus_sample = all_expr_vals
    all_expr_vals_minus_sample.remove(sample_expr_val)

    def add_noise(val):
        val = val + numpy.random.normal(0, 0.01)
        return val

    q_vals = list()
    for i in range(11):
        trial_sample_expr_val = add_noise(sample_expr_val)
        trial_expr_vals = [add_noise(x) for x in all_expr_vals_minus_sample]
        ecdf = ECDF(trial_expr_vals + [trial_sample_expr_val])
        q = ecdf(trial_sample_expr_val)
        q_vals.append(q)

    region_gene_list = ",".join(gene_list)
    # print("Q_vals: {}".format(q_vals))
    sample_q = statistics.median(q_vals)

    """
    print(
        "{}\t[{}]\t{}\t[{}]\t{}".format(
            sample_q,
            ",".join([str(x) for x in q_vals]),
            sample_expr_val,
            ",".join([str(x) for x in all_expr_vals]),
            row,
        )
    )
    """

    # print("sample_expr_val: {}".format(sample_expr_val))
    # print("mean_expr_val: {}".format(mean_expr_val))
    # print("sample_q: {}".format(sample_q))

    pseudocount = 0.01
    fold_change = math.log2(
        (sample_expr_val + pseudocount) / (mean_expr_val + pseudocount)
    )

    tab_writer.writerow(
        {
            "cohort": row["cohort"],
            "project": row["project"],
            "sample_id": row["sample_id"],
            "seqtype": row["seqtype"],
            "chrom": row["humanchr"],
            "virus": row["virus"],
            "contig": row["contig"],
            "sample_region_expr": "{:.2f}".format(sample_expr_val),
            "mean_region_expr": "{:.2f}".format(mean_expr_val),
            "log2_fold_change": "{:.2f}".format(fold_change),
            "expr_quantile": "{:.2f}".format(sample_q),
            "region_gene_list": region_gene_list,
        }
    )

    return


def get_project_bed_file_mappings(expr_bed_dir):

    files = glob.glob(expr_bed_dir + "/*sorted.gz")

    project_bed_file_mappings = dict()
    for filename in files:
        bname = os.path.basename(filename)
        if re.match("^TCGA", bname):
            tcga_tok = bname.split(".")[0].split("-")[1]
            project_bed_file_mappings[tcga_tok] = filename
        elif re.match("^GTEx", bname):
            gtex_project = bname.split(".")[1]
            project_bed_file_mappings[gtex_project] = filename

    return project_bed_file_mappings


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
