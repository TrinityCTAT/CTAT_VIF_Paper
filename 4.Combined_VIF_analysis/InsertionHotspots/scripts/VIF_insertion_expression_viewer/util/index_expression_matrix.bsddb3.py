#!/usr/bin/env python3

import sys, os, re
import bsddb3
import json
import gzip
import subprocess

usage = "usage: expr.matrix [another_expr_matrix ...] \n\n"

if len(sys.argv) < 2:
    exit(usage)

expr_matrix_files = sys.argv[1:]


def main():

    for expr_matrix_file in expr_matrix_files:

        db_fname = expr_matrix_file + ".bdb"
        dbase = bsddb3.hashopen(db_fname, "c")

        process_n_store_expr_matrix(dbase, expr_matrix_file)

        dbase.close()

    print("done")

    sys.exit(0)


def process_n_store_expr_matrix(dbase, expr_matrix_file):

    print("-processing " + expr_matrix_file)

    if re.search(".gz$", expr_matrix_file):
        fh = gzip.open(expr_matrix_file, "rt")
        linecount = subprocess.check_output(
            "gunzip -c {} | wc -l".format(expr_matrix_file), shell=True
        ).decode()
    else:
        fh = open(expr_matrix_file, "rt")
        linecount = subprocess.check_output(
            "cat {} | wc -l".format(expr_matrix_file), shell=True
        ).decode()

    linecount = int(linecount.strip())
    # print("linecount: [{}]".format(linecount))

    header = next(fh)
    header = header.rstrip()
    samples_list = header.split("\t")
    samples_list.pop(0)

    samples_list = reconstruct_sample_list_using_TCGA_names(
        expr_matrix_file, samples_list
    )

    dbase[b"sample_list"] = json.dumps(samples_list)

    print(f"sample_list for {expr_matrix_file} : {samples_list}")

    counter = 0
    for line in fh:
        counter += 1
        if counter % 100 == 0:
            sys.stderr.write("\r[{:.2f}% loaded]  ".format(counter / linecount * 100))

        line = line.rstrip()
        vals = line.split("\t")
        gene_id = vals.pop(0)
        sample_expr_list = list()
        for i, val in enumerate(vals):
            sample_name = samples_list[i]
            sample_gene_expr_struct = {
                "sample": sample_name,
                "expr": float(val),
                "rank": None,
            }

            sample_expr_list.append(sample_gene_expr_struct)

        sample_expr_list = sorted(
            sample_expr_list, key=lambda x: x["expr"], reverse=True
        )
        for i in range(len(sample_expr_list)):
            sample_expr_list[i]["rank"] = i

        for entry in sample_expr_list:
            sample_name = entry["sample"]
            sample_gene_name = sample_name + "::" + gene_id
            dbase[sample_gene_name.encode()] = json.dumps(entry)

    print("-done processing " + expr_matrix_file)

    return


def reconstruct_sample_list_using_TCGA_names(expr_matrix_filename, samples_list):

    m = re.search("TCGA-([^\.]+)\.", expr_matrix_filename)
    if not m:
        exit(
            "Sorry, cannot extract TCGA cancer type from filename: {}".format(
                expr_matrix_filename
            )
        )

    new_samples_list = list()
    TCGA_type = m.group(1)

    for sample in samples_list:
        sample = re.sub("^TCGA", TCGA_type, sample)

        if re.search("-01B$", sample) or re.search("-11A", sample):
            sample = re.sub("-[01][01][AB]$", "-NT", sample)
        elif re.search("-01A$", sample):
            sample = re.sub("-01A$", "-TP", sample)

        new_samples_list.append(sample)

    return new_samples_list


if __name__ == "__main__":
    main()
