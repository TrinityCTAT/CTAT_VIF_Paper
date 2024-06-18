#!/usr/bin/env python3

import sys, os, re

usage = f"""\n\nusage: {sys.argv[0]} ref_annot.gtf.gene_spans | \\
                        sort -k 1,1 -k2,2g -k3,3g | \\
                        bgzip -c > ref_annot.gtf.gene_spans.bed.gz && tabix ref_annot.gtf.gene_spans.bed.gz\n\n"""


if len(sys.argv) < 2:
    exit(usage)

    
gene_spans_file = sys.argv[1]


with open(gene_spans_file) as fh:
    for line in fh:
        line = line.rstrip()
        vals = line.split("\t")
        print("\t".join([*vals[1:4], vals[5] + "^" + vals[0], vals[4] ]))

