#!/usr/bin/env python3


import sys, os, re
import subprocess
import logging
import argparse

## Set up the logging  
logging.basicConfig(format='\n %(levelname)s : %(message)s', level=logging.DEBUG)
logger = logging.getLogger(__name__)


utildir = os.path.abspath(os.path.dirname(__file__))


#/home/mbrown/GitHub/VirusInsertionBenchmarking/util/Other_Pipelines/VirusFinder2/reformatVirusFinder2.py \
def reformat(virus):
    cmd = f'''
    {utildir}/reformatVirusFinder2.py \
        --predicted_insertions {virus}_results-virus-loci.txt \
        --virus {virus} \
        --sample_tag subset.insertion_seqs.fa \
        --output {virus}
    '''

    subprocess.run(cmd, shell=True)


def benchmarking(virus, insertion_truth_set):

    cmd = f'''
    {utildir}/Benchmarking.py \
            --truth_insertions {insertion_truth_set} \
            --predicted_insertions {virus}_VirusFindeer2_output.tsv \
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
    parser = argparse.ArgumentParser(description="Reformat and benchmark VirusFinder2 output for benchamrking.", 
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--Directory", type=str, required=True, help="Directory of the Virusfinder output.")
    parser.add_argument("--insertion_truth_set", type=str, required=True, help="insertion truth set file.")
    parser.add_argument("--viruses", type=list, required=False, default=["HPV45","HPV39","HPV33","HPV16","HPV31","HPV18", "HPV35"], help="Virus identified.")
    #parser.add_argument("--output", type=str, required=False, default=".", help="Directory of the ViFi output.")

    # Parse the variables given 
    args = parser.parse_args()
    Directory = args.Directory
    virus_list = args.viruses
    insertion_truth_set = args.insertion_truth_set


    # Change directory 
    os.chdir(Directory)


    # virus_list = ["HPV45","HPV39","HPV33","HPV16","HPV31","HPV18","HPV35"]

    for virus in virus_list:
        f"\trunning: {virus}"

        # Check if Directory exists, if it doesnt then make it 
        isExist = os.path.exists(virus)
        if isExist == False:
            os.makedirs(virus)

        # reformat
        reformat(virus)

        # change directory
        os.chdir(f'{virus}')
        # benchamrking
        benchmarking(virus, insertion_truth_set)
        TP_FP()

        # change directory again 
        os.chdir(Directory)


if __name__ == "__main__":
    main()
