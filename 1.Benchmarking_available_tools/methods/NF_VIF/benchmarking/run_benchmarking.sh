#!/bin/bash
set -e


./scripts/benchmark_wrapper.py \
    --output_directory `pwd`/data \
    --viruses HPV45,HPV39,HPV33,HPV35,HPV16,HPV31,HPV18 \
    --output `pwd`/data/combined_vf_VIF_output.tsv
    

./scripts/Benchmarking.py \
        --truth_insertions `pwd`/insertion_truth_set.tsv \
        --predicted_insertions `pwd`/data//combined_vf_VIF_output.tsv \
        --output `pwd`/benchmarking_out


./scripts/analyze_VIF_TP_FP_FNs.Rscript \
    --dat `pwd`/benchmarking_out/Benchmarking_Output.txt



