#!/bin/bash




set -e



./Benchmarking/run_benchmarking_VB.py \
--Directory `pwd` \
--insertion_truth_set /home/mbrown/CTAT_VIF/Simulation/CreateSimulations_50/insertion_truth_set.tsv
