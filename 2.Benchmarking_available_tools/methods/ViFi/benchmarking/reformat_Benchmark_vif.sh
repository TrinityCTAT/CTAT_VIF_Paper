#!/bin/bash

set -ex


# Reformat
./scripts/reformatViFi.py \
        --inputs ./data/Test.clusters.txt \
        --sampleID subset.insertion_seqs.fa \
        --output `pwd`/benchmarking_out


./scripts/Benchmarking.py \
        --truth_insertions `pwd`/insertion_truth_set.tsv \
        --predicted_insertions `pwd`/benchmarking_out/ViFi_Output_reformated.txt \
        --output `pwd`/benchmarking_out

# STATS
./scripts/analyze_VIF_TP_FP_FNs.Rscript \
        --dat `pwd`/benchmarking_out/Benchmarking_Output.txt



