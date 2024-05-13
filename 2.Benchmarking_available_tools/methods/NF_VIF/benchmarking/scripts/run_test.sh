#!/bin/bash

set -e




./run_wrapper.py \
	--fastq1 /home/mbrown/CTAT_VIF/Simulation/CreateSimulations_50/subset.insertion_seqs_QSadjust_left.fq \
	--fastq2 /home/mbrown/CTAT_VIF/Simulation/CreateSimulations_50/subset.insertion_seqs_QSadjust_right.fq \
	--viruses HPV45,HPV39,HPV33,HPV35,HPV16,HPV31,HPV18 \
	--output /mnt/disks/ses-extra/mbrown/CTAT_VIF/NF_VIF/OUTPUT_Simulation


cd ../OUTPUT_Simulation/HPV45
chmod 744 run_command.sh 
./run_command.sh


cd ../HPV39
chmod 744 run_command.sh 
./run_command.sh

cd ../HPV33
chmod 744 run_command.sh 
./run_command.sh

cd ../HPV35
chmod 744 run_command.sh 
./run_command.sh

cd ../HPV16
chmod 744 run_command.sh 
./run_command.sh

cd ../HPV31
chmod 744 run_command.sh 
./run_command.sh

cd ../HPV18
chmod 744 run_command.sh 
./run_command.sh