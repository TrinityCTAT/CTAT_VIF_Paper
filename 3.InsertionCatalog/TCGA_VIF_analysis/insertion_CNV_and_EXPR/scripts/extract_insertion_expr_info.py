#!/usr/bin/env python3

import sys, os, re
import argparse
import json
import bsddb3
from collections import defaultdict
import csv
import subprocess


def main():

    parser = argparse.ArgumentParser(
        description="extract expr view",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "--insertions_tsv",
        type=str,
        required=True,
        help="insertions input file (requires neighbors set)",
    )

    parser.add_argument(
        "--expr_matrix_bdbs",
        type=str,
        required=True,
        help="gene expr matrix stored in pickle format",
        nargs="+",
    )

    parser.add_argument(
        "--output_filename",
        type=str,
        required=True,
        help="name of output file (can include pdf or svg extension",
    )

    args = parser.parse_args()

    expr_matrix_bdbs = args.expr_matrix_bdbs
    vif_insertions_tsv = args.insertions_tsv
    output_filename = args.output_filename

    # define TCGA class to bdb file.
    def TCGA_to_bdb(bdbs_list):
        tcga_class_to_bdb_filename = dict()
        for bdb in bdbs_list:
            bname = os.path.basename(bdb)
            m = re.search("^TCGA-([^\\.]+)", bname)
            if not m:
                raise RuntimeError("cannot decipher TCGA class from {}".format(bname))
            tcga_class = m.group(1)
            tcga_class_to_bdb_filename[tcga_class] = bdb

        return tcga_class_to_bdb_filename

    tcga_class_to_bdb = TCGA_to_bdb(expr_matrix_bdbs)
    print(tcga_class_to_bdb)

    ofh = open(output_filename, "wt")
    fh = open(vif_insertions_tsv, "rt")

    csv_reader = csv.DictReader(fh, delimiter="\t")
    outfields = ["TCGA", "sample", "virus", "insertion", "gene", "expr", "rank"]
    csv_writer = csv.DictWriter(ofh, outfields, delimiter="\t")
    csv_writer.writeheader()

    row_counter = 0

    for row in csv_reader:

        row_counter += 1
        if row_counter % 10 == 0:
            sys.stderr.write("\r{}  ".format(row_counter))

        sample = row["sample"]
        tcga = sample.split("-")[0]
        virus = row["virus"]
        insertion = row["contig"]
        genes = (
            row["left_genes"].split(";")
            + row["insertion_genes"].split(";")
            + row["right_genes"].split(";")
        )

        bdb = tcga_class_to_bdb[tcga]

        for gene in genes:
            expr, rank = get_gene_expr_info(sample, gene, bdb)
            # print(f"{sample}::{gene} = {expr} {rank}")
            if expr is not None:
                csv_writer.writerow(
                    {
                        "TCGA": tcga,
                        "sample": sample,
                        "virus": virus,
                        "insertion": insertion,
                        "gene": gene,
                        "expr": f"{expr:0.3f}",
                        "rank": rank,
                    }
                )

    print("\n\nDone.\n\n")

    sys.exit(0)


def get_gene_expr_info(sample, target_gene, expr_matrix_bdb):

    dbase = bsddb3.hashopen(expr_matrix_bdb, "r")
    all_samples = set(json.loads(dbase[b"sample_list"]))
    all_genes = set(json.loads(dbase[b"gene_list"]))

    if sample not in all_samples:
        return (None, None)

    symbol_to_gene_id = defaultdict(str)
    for gene in all_genes:
        symbol = gene.split("^")[0]
        symbol_to_gene_id[symbol] += gene

    if target_gene in symbol_to_gene_id:
        full_gene_id = symbol_to_gene_id[target_gene]
        if full_gene_id in all_genes:
            sample_gene_key = sample + "::" + full_gene_id
            expr_info = json.loads(dbase[sample_gene_key.encode()])
            expr = expr_info["expr"]
            rank = expr_info["rank"]

            return (expr, rank)

    # not found
    return (None, None)


if __name__ == "__main__":
    main()
