#!/usr/bin/env python3


import sys, os, re
import subprocess
import logging
import argparse

## Set up the logging  
logging.basicConfig(format='\n %(levelname)s : %(message)s', level=logging.DEBUG)
logger = logging.getLogger(__name__)


utildir = os.path.abspath(os.path.dirname(__file__))


def reformat(virus):
	cmd = f'''
	{utildir}/reformatingVIRUSBreakend.py \
		--predicted_insertions cromwell-executions/VIRUSBreakend/*/call-RunVIRUSBreakend/execution/sample.virusbreakend.vcf \
		--virus {virus} \
		--sample_tag subset.insertion_seqs.fa \
		--output "."
	'''

	# print(cmd)

	subprocess.run(cmd, shell=True)


def benchmarking(virus, insertion_truth_set):

	cmd = f'''
	{utildir}/Benchmarking.py \
	        --truth_insertions {insertion_truth_set} \
	        --predicted_insertions {virus}_VIRUSBreakend_output.tsv \
			--Viruses {virus}
	'''
	subprocess.run(cmd, shell=True)

def TP_FP():
	cmd = f'''
	{utildir}/analyze_VIF_TP_FP_FNs.Rscript \
		--dat Benchmarking_Output.txt
	'''
	subprocess.run(cmd, shell=True)




def main():

	####################
    # Parse the use supplied information 
    ####################
    # Set the variables to supply 
    parser = argparse.ArgumentParser(description="Reformat and benchmark Virusbreakend output for benchamrking.", 
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--Directory", type=str, required=True, help="Directory of the Virusbreakend output.")
    parser.add_argument("--insertion_truth_set", type=str, required=True, help="insertion truth set file.")
    parser.add_argument("--viruses", type=list, required=False, default=["HPV45","HPV39","HPV33","HPV16","HPV31","HPV18", "HPV35"], help="Virus identified.")
    #parser.add_argument("--output", type=str, required=False, default=".", help="Directory of the ViFi output.")

    # Parse the variables given 
    args = parser.parse_args()
    Directory = args.Directory
    virus_list = args.viruses
    insertion_truth_set = args.insertion_truth_set



    os.chdir(Directory)

    # virus_list = ["HPV45","HPV39","HPV33","HPV16","HPV31","HPV18", "HPV35"]

    for virus in virus_list:
    	f"\trunning: {virus}"
    	# change directory
    	os.chdir(f'{virus}')

    	# reformat
    	reformat(virus)
    	# benchamrking
    	benchmarking(virus, insertion_truth_set)
    	TP_FP()

    	# change directory again 
    	os.chdir(Directory)


if __name__ == "__main__":
    main()
