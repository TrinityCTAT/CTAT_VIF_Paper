#!/usr/bin/env python3

import sys, os, re
from collections import defaultdict
import csv
import argparse
import pandas as pd


def main():

    parser = argparse.ArgumentParser(
        description="assign neighboring genes to hotspots",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "--hotspots",
        type=str,
        required=True,
        help="filename hotspot_sample_data.threshold=NUM.tsv",
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
    parser.add_argument(
        "--no_gene_decoration",
        action="store_true",
        default=False,
        help="do not highlight any gene insertion",
    )

    parser.add_argument(
        "--only_top_insertion_per_sample_hotspot",
        action="store_true",
        default=False,
        help="restrict output to topmost supported insertion per sample hotspot",
    )

    parser.add_argument("--output", type=str, required=True, help="output tsv filename")

    args = parser.parse_args()

    hotspots_file = args.hotspots
    ref_gene_spans_file = args.ref_gene_spans
    num_genes_include = args.num_genes_include
    no_gene_decoration_flag = args.no_gene_decoration
    output_filename = args.output
    only_top_insertion_per_sample_hotspot = args.only_top_insertion_per_sample_hotspot

    hotspots_df = pd.read_csv(hotspots_file, sep="\t")

    # want the single top insertion site per sample/hotspot combo.

    if only_top_insertion_per_sample_hotspot:

        def get_top_insertion(df):
            df = df.sort_values(["total"], ascending=False)
            return df.head(1)

        hotspots_df = hotspots_df.groupby(["sample", "hotspot"]).apply(
            get_top_insertion
        )
        hotspots_df = hotspots_df.reset_index(drop=True)

    chr_to_ordered_gene_list = parse_gene_spans(ref_gene_spans_file)

    def annotate_neighboring_genes_for_row(row):
        chr = row["human_chrom"]
        coord = row["human_coord"]

        left_genes, insert_genes, right_genes = annotate_neighboring_genes(
            chr,
            coord,
            chr_to_ordered_gene_list[chr],
            num_genes_include,
            no_gene_decoration_flag,
        )

        row["left_genes"] = ";".join(left_genes)
        row["insert_genes"] = ";".join(insert_genes)
        row["right_genes"] = ";".join(right_genes)

        return row

    hotspots_df = hotspots_df.apply(annotate_neighboring_genes_for_row, axis=1)

    hotspots_df["hotspot_coord"] = hotspots_df["Start_b"]

    hotspots_df = hotspots_df.sort_values(["hotspot"])

    """
    hotspots_df = hotspots_df[
        [
            "sample",
            "cohort",
            "project",
            "contig",
            "human_chrom",
            "human_coord",
            "virus_genome",
            "virus_coord",
            "virus_brkend_grp",
            "is_primary",
            "splice_type",
            "total",
            "total_rpm",
            "hotspot",
            "hotspot_coord",
            "hotspot_sample_counts",
            "left_genes",
            "insert_genes",
            "right_genes",
        ]
    ]
    """

    hotspots_df.to_csv(output_filename, index=False, sep="\t")

    sys.exit(0)


def annotate_neighboring_genes(
    chrom, coord, ordered_gene_list, num_genes_include, no_gene_decoration_flag
):

    left_genes = list()
    insert_genes = list()
    right_genes = list()

    coord = int(coord)

    for gene in ordered_gene_list:
        if gene["rend"] < coord:
            left_genes.append(gene["gene_sym"])
        elif gene["lend"] > coord:
            right_genes.append(gene["gene_sym"])
        elif gene["lend"] <= coord and gene["rend"] >= coord:
            insert_genes.append(gene["gene_sym"])

    # want last few left genes
    left_genes = left_genes[-1 * num_genes_include :]

    left_genes.extend(insert_genes)

    # want first few right genes
    right_genes = right_genes[0 : min(len(right_genes), num_genes_include)]

    return left_genes, insert_genes, right_genes


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
