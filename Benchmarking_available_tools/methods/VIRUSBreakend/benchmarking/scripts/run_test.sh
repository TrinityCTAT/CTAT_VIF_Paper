#!/bin/bash

set -e




# ./run_wrapper.py \
# 	--fastq1 /home/mbrown/CTAT_VIF/Simulation/CreateSimulations_2/subset.insertion_seqs_QSadjust_left.fq \
# 	--fastq2 /home/mbrown/CTAT_VIF/Simulation/CreateSimulations_2/subset.insertion_seqs_QSadjust_right.fq \
# 	--viruses HPV45,HPV39,HPV33,HPV35,HPV16,HPV31,HPV18 \
# 	--output /mnt/disks/ses-extra/mbrown/CTAT_VIF/VIRUSBreakend



./run_wrapper.py \
	--fastq1 /home/mbrown/CTAT_VIF/Simulation/CreateSimulations_50/subset.insertion_seqs_QSadjust_50ins_7sam_left.fq \
	--fastq2 /home/mbrown/CTAT_VIF/Simulation/CreateSimulations_50/subset.insertion_seqs_QSadjust_50ins_7sam_right.fq \
	--viruses HPV45,HPV39,HPV33,HPV35,HPV16,HPV31,HPV18 \
	--output /mnt/disks/ses-extra/mbrown/CTAT_VIF/VIRUSBreakend



