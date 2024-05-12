#!/bin/bash

set -e




./run_wrapper.py \
	--fastq1 /home/mbrown/CTAT_VIF/Simulation/CreateSimulations_50/subset.insertion_seqs_QSadjust_50ins_7sam_left.fq \
	--fastq2 /home/mbrown/CTAT_VIF/Simulation/CreateSimulations_50/subset.insertion_seqs_QSadjust_50ins_7sam_right.fq \
	--viruses HPV45,HPV39,HPV33,HPV35,HPV16,HPV31,HPV18 \
	--output `pwd`/TEST

cd TEST

#~~~~~~~ Integrations ~~~~~~~~~~
# mkdir Output


# Run each 
for element in HPV45 HPV39 HPV33 HPV35 HPV16 HPV31 HPV18
do
	cd $element

	FILE=cromwell-executions
	if test -f "$FILE"; then
	    sudo rm -r cromwell-*
	fi

	chmod 744 run_command.sh
	./run_command.sh

	cd ..

done