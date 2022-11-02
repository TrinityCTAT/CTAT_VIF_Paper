#!/usr/bin/env python3

import sys, os, re
from collections import defaultdict
import csv
import argparse
import pandas as pd


def main():

    parser = argparse.ArgumentParser(
        description="assign neighboring genes to insertions",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "--insertions", type=str, required=True, help="filename insertions.tsv",
    )
    parser.add_argument(
        "--ref_gene_spans", type=str, required=True, help="ref annot gene spans file"
    )
    parser.add_argument(
        "--num_genes_include",
        type=int,
        default=3,
        help="number of genes to include on both sides of the breakpoint",
    )

    parser.add_argument("--output", type=str, required=True, help="output tsv filename")

    args = parser.parse_args()

    insertions_file = args.insertions
    ref_gene_spans_file = args.ref_gene_spans
    num_genes_include = args.num_genes_include
    output_filename = args.output

    insertions_df = pd.read_csv(insertions_file, sep="\t")

    chr_to_ordered_gene_list = parse_gene_spans(ref_gene_spans_file)

    def annotate_neighboring_genes_for_row(row):
        chr = row["human_chrom"]
        coord = row["human_coord"]

        left_genes, insertion_genes, right_genes = annotate_neighboring_genes(
            chr, coord, chr_to_ordered_gene_list[chr], num_genes_include,
        )

        row["left_genes"] = ";".join(left_genes)
        row["insertion_genes"] = ";".join(insertion_genes)
        row["right_genes"] = ";".join(right_genes)

        return row

    insertions_df = insertions_df.apply(annotate_neighboring_genes_for_row, axis=1)

    insertions_df.to_csv(output_filename, index=False, sep="\t")

    sys.exit(0)


def annotate_neighboring_genes(chrom, coord, ordered_gene_list, num_genes_include):

    left_genes = list()
    right_genes = list()

    insertion_genes = list()

    coord = int(coord)

    for gene in ordered_gene_list:
        if gene["rend"] < coord:
            left_genes.append(gene["gene_sym"])
        elif gene["lend"] > coord:
            right_genes.append(gene["gene_sym"])
        elif gene["lend"] <= coord and gene["rend"] >= coord:
            insertion_genes.append(gene["gene_sym"])

    # want last few left genes
    left_genes = left_genes[-1 * num_genes_include :]

    # want first few right genes
    right_genes = right_genes[0 : min(len(right_genes), num_genes_include)]

    return left_genes, insertion_genes, right_genes


def parse_gene_spans(ref_gene_spans_file):

    chrom_to_gene_list = defaultdict(list)

    with open(ref_gene_spans_file) as fh:
        for line in fh:
            line = line.rstrip()
            ensg_id, chrom, lend, rend, orient, gene_sym, gene_type = line.split("\t")

            feature = {
                "chrom": chrom,
                "lend": int(lend),
                "rend": int(rend),
                "midpt": int((int(lend) + int(rend)) / 2),
                "orient": orient,
                "gene_sym": gene_sym,
                "gene_type": gene_type,
            }

            chrom_to_gene_list[chrom].append(feature)

    for chrom, gene_list in chrom_to_gene_list.items():
        gene_list = sorted(gene_list, key=lambda x: x["midpt"])
        chrom_to_gene_list[chrom] = gene_list

    return chrom_to_gene_list


if __name__ == "__main__":
    main()
