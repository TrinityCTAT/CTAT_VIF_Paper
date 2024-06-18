#!/usr/bin/env python

import sys, os, re

usage = "\n\n\tusage:  TCGA-...htseq_fpkm-up.tsv  ref_annot.gtf.gene_spans >  updated.matrix\n\n"

if len(sys.argv) < 3:
    exit(usage)

    
fpkm_uq_matrix = sys.argv[1]
ref_annot_gtf_gene_spans = sys.argv[2]

gene_sym_map = dict()

with open(ref_annot_gtf_gene_spans) as fh:
    for line in fh:
        line = line.rstrip()
        vals = line.split("\t")
        ensg_id = vals[0]
        gene_sym = vals[5]

        gene_sym_map[ensg_id] = gene_sym

with open(fpkm_uq_matrix) as fh:
    header = next(fh)
    header = header.rstrip()
    print(header)
    for line in fh:
        line = line.rstrip()
        vals = line.split("\t")
        ensg_id = vals[0]
        if ensg_id in gene_sym_map:
            gene_sym = gene_sym_map[ensg_id]
            vals[0] = gene_sym + "^" + ensg_id
        print("\t".join(vals))


