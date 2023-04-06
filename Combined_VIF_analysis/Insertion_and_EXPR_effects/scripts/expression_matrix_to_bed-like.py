#!/usr/bin/env python3

import sys, os, re
import gzip
import subprocess


def main():

    usage = "\n\tusage: {} expression_matrix.gz gene_spans_file\n\n".format(sys.argv[0])

    if len(sys.argv) < 3:
        sys.exit(usage)

    expr_matrix_file = sys.argv[1]
    gene_spans_file = sys.argv[2]

    TCGA_token = os.path.basename(expr_matrix_file).split(".")[0].split("-")[1]

    gene_spans_info = parse_gene_spans(gene_spans_file)

    expr_matrix = gzip.open(expr_matrix_file, "rt")
    header = next(expr_matrix)
    header = header.rstrip()
    samples = header.split("\t")[1:]

    samples = [revise_sample_name(x, TCGA_token) for x in samples]

    with open(expr_matrix_file + ".sample_names", "wt") as ofh:
        print(",".join(samples), file=ofh)

    missing_counter = 0

    gene_expr_gene_span_list_filename = expr_matrix_file + ".gene_expr_gene_span_list"

    with open(gene_expr_gene_span_list_filename, "wt") as ofh:

        for row in expr_matrix:
            row = row.rstrip()
            row = row.split("\t")
            gene = row.pop(0)
            vals = ["{:.2f}".format(float(x)) for x in row]

            if gene in gene_spans_info:
                gene_span_info = gene_spans_info[gene]

                print("\t".join([gene_span_info, gene, ",".join(vals)]), file=ofh)
            else:
                print("-missing span info for gene: " + gene, file=sys.stderr)
                missing_counter += 1

    print("\nMissing gene spans for {} entries.".format(missing_counter))

    cmd = "sort -k1,1 -k2,2g -k3,3g {} > {}.sorted".format(
        gene_expr_gene_span_list_filename, gene_expr_gene_span_list_filename
    )
    subprocess.check_call(cmd, shell=True)

    cmd = "bgzip -f {}.sorted".format(gene_expr_gene_span_list_filename)
    subprocess.check_call(cmd, shell=True)

    cmd = "tabix -f -p bed {}.sorted.gz".format(gene_expr_gene_span_list_filename)
    subprocess.check_call(cmd, shell=True)

    sys.exit(0)


def revise_sample_name(sample_name, TCGA_token):

    sample_name = sample_name.replace("TCGA", TCGA_token)

    if re.search("-01B$", sample_name) or re.search("-11A", sample_name):
        sample_name = re.sub("-[01][01][AB]$", "-NT", sample_name)
    elif re.search("-01A$", sample_name):
        sample_name = re.sub("-01A$", "-TP", sample_name)

    return sample_name


def parse_gene_spans(gene_spans_file):

    gene_spans_info = dict()

    with open(gene_spans_file, "rt") as fh:
        for line in fh:
            line = line.rstrip()
            vals = line.split("\t")
            ensmbl_id = vals[0]
            chrom = vals[1]
            lend = vals[2]
            rend = vals[3]
            gene_symbol = vals[5]

            gene_id = gene_symbol + "^" + ensmbl_id

            gene_spans_info[gene_id] = "\t".join([chrom, lend, rend])

    return gene_spans_info


if __name__ == "__main__":
    main()
