#!/bin/bash
set -e








./run_benchmarking_VF2.py \
    --Directory /mnt/disks/ses-extra/mbrown/CTAT_VIF/VirusFinder2/Output \
    --insertion_truth_set /home/mbrown/CTAT_VIF/Simulation/CreateSimulations_50/insertion_truth_set.tsv
