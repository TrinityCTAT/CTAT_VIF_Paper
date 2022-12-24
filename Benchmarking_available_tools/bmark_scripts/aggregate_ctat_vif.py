#!/usr/bin/env python

import sys, os, re

usage = f"usage: {sys.argv[0]} fileA [ fileB ...]\n\n"
if len(sys.argv) < 2:
    exit(usage)


files = sys.argv[1:]
got_header = False
for file in files:
    with open(file) as fh:
        header = next(fh)
        if not got_header:
            got_header = True
            sys.stdout.write("sample_name\t" + header)

        sample_name = os.path.basename(file)
        sample_name = sample_name.replace(".vif.refined.tsv", "")
        for line in fh:
            sys.stdout.write(f"{sample_name}\t" + line)

sys.exit(0)
