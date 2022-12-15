#!/usr/bin/env python3


import sys, os, re
import subprocess
import logging
import argparse
import glob
import pandas as pd

## Set up the logging  
logging.basicConfig(format='\n %(levelname)s : %(message)s', level=logging.DEBUG)
logger = logging.getLogger(__name__)

utildir = os.path.abspath(os.path.dirname(__file__))


#/home/mbrown/GitHub/VirusInsertionBenchmarking/util/Other_Pipelines/VirusFinder2/reformatVirusFinder2.py \
def reformat(virus):
    # Get the current working directory
    cwd = os.getcwd()


    logger.info(f"\t Current Virus: {virus}")
    # output directory path 
    output_dir = os.path.join(cwd,virus)
    # input file of interest 
    ## Search in the directory 
    files_list = glob.glob(f"{output_dir}/cromwell-executions/*/*/*/*/*")
    infile  = [i for i in files_list if i.split("/")[-1] == "final_hits.txt"][0]

    cmd = f'''
    {utildir}/reformatBatVI.py \
        --predicted_insertions {infile} \
        --virus {virus} \
        --output {output_dir} \
        --sample_tag subset.insertion_seqs.fa
        '''
    subprocess.run(cmd, shell=True)

def combine(virus_list):
    # Get the current working directory
    cwd = os.getcwd()
    df = pd.DataFrame()
    for virus in virus_list:
        input_file = f"{virus}_BatVI_output.tsv"

        Batvi_df =  pd.read_csv(input_file,sep = "\t")

        df = pd.concat([df, Batvi_df])

    df = df.reset_index()
    df.to_csv(f"{cwd}/combined_BatVI_output.tsv", sep = "\t", index=False)


def benchmarking(insertion_truth_set, Directory):

    cmd = f'''
    {utildir}/Benchmarking.py \
        --truth_insertions {insertion_truth_set} \
        --predicted_insertions combined_BatVI_output.tsv \
        --output {Directory}
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
    parser.add_argument("--Directory", type=str, required=True, help="Directory of the BATVI output.")
    parser.add_argument("--insertion_truth_set", type=str, required=True, help="insertion truth set file.")
    parser.add_argument("--viruses", type=list, required=False, default=["HPV45","HPV39","HPV33","HPV16","HPV31","HPV18", "HPV35"], help="Virus identified.")
    #parser.add_argument("--output", type=str, required=False, default=".", help="Directory of the ViFi output.")

    # Parse the variables given 
    args = parser.parse_args()
    Directory = args.Directory
    virus_list = args.viruses
    insertion_truth_set = args.insertion_truth_set


    os.chdir(Directory)

    for virus in virus_list:
        f"\trunning: {virus}"
        # change directory
        # os.chdir(f'{virus}')
        # reformat
        reformat(virus)

        
    combine(virus_list)
        # benchamrking
    benchmarking(insertion_truth_set, Directory)
    TP_FP()



if __name__ == "__main__":
    main()
