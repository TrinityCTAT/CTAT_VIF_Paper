#!/bin/bash
set -e


./benchamrk_wrapper.py \
    --output_directory /mnt/disks/ses-extra/mbrown/CTAT_VIF/NF_VIF/OUTPUT_Simulation \
    --viruses HPV45,HPV39,HPV33,HPV35,HPV16,HPV31,HPV18 \
    --output /mnt/disks/ses-extra/mbrown/CTAT_VIF/NF_VIF/Benchmarking/combined_vf_VIF_output.tsv




/home/mbrown/GitHub/VirusInsertionBenchmarking_old/VirusInsertionBenchmarking/util/Other_Pipelines/Benchmarking.py \
        --truth_insertions /home/mbrown/CTAT_VIF/Simulation/CreateSimulations_50/insertion_truth_set.tsv \
        --predicted_insertions /mnt/disks/ses-extra/mbrown/CTAT_VIF/NF_VIF/Benchmarking/combined_vf_VIF_output.tsv \
        --output /mnt/disks/ses-extra/mbrown/CTAT_VIF/NF_VIF/Benchmarking





/home/mbrown/GitHub/VirusInsertionBenchmarking_old/VirusInsertionBenchmarking/util/analyze_VIF_TP_FP_FNs.Rscript \
    --dat /mnt/disks/ses-extra/mbrown/CTAT_VIF/NF_VIF/Benchmarking/Benchmarking_Output.txt
