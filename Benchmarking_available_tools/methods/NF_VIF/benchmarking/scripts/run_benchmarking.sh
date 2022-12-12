#!/bin/bash
set -e


# ./Reformating_vf_VIF.py \
#     --output_directory /home/mbrown/CTAT_VIF/Other_Pipelines/NF_VIF/output/hpv_mapping/blat \
#     --sample_tag subset.insertion_seqs.fa \
#     --output /home/mbrown/CTAT_VIF/Other_Pipelines/NF_VIF/Benchmarking/vf_VIF_output.tsv

./benchamrk_wrapper.py \
    --output_directory /home/mbrown/CTAT_VIF/Other_Pipelines/NF_VIF/OUTPUT_Simulation \
    --viruses HPV45,HPV39,HPV33,HPV35,HPV16,HPV31,HPV18 \
    --output /home/mbrown/CTAT_VIF/Other_Pipelines/NF_VIF/OUTPUT_Simulation/combined_vf_VIF_output.tsv
    
# ./Reformating_vf_VIF.py \
#     --output_directory /home/mbrown/CTAT_VIF/Other_Pipelines/NF_VIF/OUTPUT_Simulation \
#     --sample_tag subset.insertion_seqs.fa \
#     --output /home/mbrown/CTAT_VIF/Other_Pipelines/NF_VIF/Benchmarking/combined_vf_VIF_output.tsv



/home/mbrown/GitHub/VirusInsertionBenchmarking_old/VirusInsertionBenchmarking/util/Other_Pipelines/Benchmarking.py \
        --truth_insertions /home/mbrown/CTAT_VIF/Simulation/CreateSimulations_50/insertion_truth_set.tsv \
        --predicted_insertions /home/mbrown/CTAT_VIF/Other_Pipelines/NF_VIF/OUTPUT_Simulation/combined_vf_VIF_output.tsv \
        --output /home/mbrown/CTAT_VIF/Other_Pipelines/NF_VIF/Benchmarking



/home/mbrown/GitHub/VirusInsertionBenchmarking_old/VirusInsertionBenchmarking/util/analyze_VIF_TP_FP_FNs.Rscript \
    --dat /home/mbrown/CTAT_VIF/Other_Pipelines/NF_VIF/Benchmarking/Benchmarking_Output.txt