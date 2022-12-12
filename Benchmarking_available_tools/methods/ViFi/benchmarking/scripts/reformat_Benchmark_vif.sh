#!/bin/bash

set -e


# Reformat
/home/mbrown/GitHub/VirusInsertionBenchmarking/util/Other_Pipelines/ViFi/reformatViFi.py \
        --inputs ./Test/*.clusters.txt \
        --sampleID subset.insertion_seqs.fa \
        --output `pwd`/benchmarking



# Benchmark
#./../Benchmarking.py
#/home/mbrown/GitHub/VirusInsertionBenchmarking/util/Other_Pipelines/Benchmarking.py \
#	--truth_insertions /home/mbrown/CTAT_VIF/Simulation/CreateSimulations/truthset.tsv \
#        --predicted_insertions /home/mbrown/CTAT_VIF/Other_Pipelines/ViFi/ViFi_Output_reformated.txt \
#        --output /home/mbrown/CTAT_VIF/Other_Pipelines/ViFi

# /home/mbrown/GitHub/VirusInsertionBenchmarking/util/Other_Pipelines/Benchmarking.py \
/home/mbrown/GitHub/VirusInsertionBenchmarking_old/VirusInsertionBenchmarking/util/Other_Pipelines/Benchmarking.py \
        --truth_insertions /home/mbrown/CTAT_VIF/Simulation/CreateSimulations_50/insertion_truth_set.tsv \
        --predicted_insertions `pwd`/benchmarking/ViFi_Output_reformated.txt \
        --output `pwd`/benchmarking

# STATS
# /home/mbrown/GitHub/VirusInsertionBenchmarking/util/analyze_VIF_TP_FP_FNs.Rscript \
/home/mbrown/GitHub/VirusInsertionBenchmarking_old/VirusInsertionBenchmarking/util/analyze_VIF_TP_FP_FNs.Rscript \
        --dat `pwd`/benchmarking/Benchmarking_Output.txt
