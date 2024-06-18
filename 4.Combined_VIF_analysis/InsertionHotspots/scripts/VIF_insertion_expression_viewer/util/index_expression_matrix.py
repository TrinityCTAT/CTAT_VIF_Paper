#!/usr/bin/env python3

import sys, os, re
import pickle
import json
from collections import defaultdict

usage = "usage: expr.matrix\n\n"

if len(sys.argv) < 2:
    exit(usage)
    
expr_matrix_file = sys.argv[1]

#gene_to_expr_matrix = dict()
sample_to_gene_expr = defaultdict(dict)


with open(expr_matrix_file) as fh:
    header = next(fh)
    header = header.rstrip()
    samples_list = header.split("\t")
    samples_list.pop(0)
    for line in fh:
        line = line.rstrip()
        vals = line.split("\t")
        gene_id = vals.pop(0)
        sample_expr_list = list()
        for i, val in enumerate(vals):
            sample_name = samples_list[i]
            sample_gene_expr_struct = {'sample':sample_name,
                                     'expr':float(val),
                                     'rank':None}
            
            sample_expr_list.append(sample_gene_expr_struct)
            sample_to_gene_expr[sample_name][gene_id] = sample_gene_expr_struct

        sample_expr_list = sorted(sample_expr_list, key=lambda x: x['expr'], reverse=True)
        for i in range(len(sample_expr_list)):
            sample_expr_list[i]['rank'] = i
            

        
        
pickle_file = expr_matrix_file + ".pickle"
with open(pickle_file, 'wb') as ofh:
    pickle.dump(sample_to_gene_expr, ofh)

print("done")

sys.exit(0)

