#!/bin/bash

set -e


export PATH=/home/mbrown/GitHub/wgsim-trans/:$PATH

VirusDB="/home/mbrown/CTAT_VIF/Simulation/CreateSimulations_50/HPV_subset.fasta"
CTAT_db="/home/bhaas/CTAT_GENOME_LIBS/GRCh38_gencode_v22_CTAT_lib_Mar012021.plug-n-play/ctat_genome_lib_build_dir"

/home/mbrown/GitHub/VirusInsertionBenchmarking/util/simulate_virus_insertions.py \
	--virus_db $VirusDB \
	--genome_lib_dir $CTAT_db \
	--ins_per_virus 50 \
	--out_prefix "subset" \
	--sim_read_length 70 2>&1 | tee $file.log




# Create the truth set
/home/mbrown/GitHub/VirusInsertionBenchmarking/util/make_truth_list.py subset > insertion_truth_set.tsv
