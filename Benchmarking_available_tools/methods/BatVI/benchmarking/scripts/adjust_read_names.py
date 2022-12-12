#!/usr/bin/env python3




import pandas as pd
import pysam 
import os, sys, re
import logging
import sys, time
import itertools
import subprocess
import argparse
logging.basicConfig(format='\n %(levelname)s : %(message)s', level=logging.DEBUG)
logger = logging.getLogger(__name__)

exec(open('FastFileClass.py').read())


def adjustNames(dic1, dic2):
    '''
    Function to reasign numeric IDs to read names in simulated data.
        Returns the relationship table to connect the numeric id to the true id 
    '''

    ID_sample = {}
    new_dic_1 = {}
    new_dic_2 = {}
    if len(dic1) != len(dic2):
        print("Error")

    # print(dic1)

    IDS = list(range(0,len(dic1)))
    keys_list = [i for i in dic1]

    if len(keys_list[0].split("/")) == 2:

        for i,j in zip(keys_list, IDS):
            # Hold Neew IDS in relasion table 
            read_id = i.split("/")[0]
            ID_sample[j] = read_id
            # new samples
            new_dic_1[j] = dic1[f"{read_id}/1"]
            new_dic_2[j] = dic2[f"{read_id}/2"]
        return new_dic_1, new_dic_2, ID_sample
    
    else:
        # for i,j in zip(keys_list, IDS):
        for i,j,z in zip(dic1, dic2, IDS):    
            # remove the tag on the end
            read_id1 = i.split(" ")[0]
            read_id2 = j.split(" ")[0]

            #check
            if read_id1 != read_id2:
                print("error")

            # Hold Neew IDS in relasion table 
            ID_sample[z] = read_id1
            # new samples
            new_dic_1[z] = dic1[i]
            new_dic_2[z] = dic2[j]
        return new_dic_1, new_dic_2, ID_sample



def rewriteFiles(left, right, output):
    #~~~~~~~~~~~~~~~~
    # Left
    #~~~~~~~~~~~~~~~~
    # Read in the fastQ
    inputfile_left = fastFile(input_file = left)
    inputfile_left = inputfile_left.fqReader()
    #~~~~~~~~~~~~~~~~
    # Right 
    #~~~~~~~~~~~~~~~~
    # Read in the fastQ
    inputfile_right = fastFile(input_file = right)
    inputfile_right = inputfile_right.fqReader()


    #~~~~~~~~~~~
    # Now replace the old sequences 
    #~~~~~~~~~~~
    new_sequences_1, new_sequences_2, relation = adjustNames(inputfile_left.sequences, inputfile_right.sequences)
    inputfile_left.sequences = new_sequences_1
    inputfile_right.sequences = new_sequences_2

    # CHECK 
    if not len(inputfile_left.sequences) == len(inputfile_right.sequences):
        m = f"\t NOT All reads found! {len(inputfile_left.sequences)} out of {len(inputfile_right.sequences)} Found."
        logger.warning(m)


    ## save the relasionship table so can access the read names later 
    relation = pd.DataFrame({"reads":relation})
    relation.to_csv(f"{output}/relationship_table.txt",sep = "\t")

    # now rewrite the new fastq file 
    new_File = f"{output}/noCHR_left.fq"
    inputfile_left.fqWriter( output = new_File)
    # now rewrite the new fastq file 
    new_File = f"{output}/noCHR_right.fq"
    inputfile_right.fqWriter( output = new_File)


def main():
    #add options to inputs
    parser = argparse.ArgumentParser(   formatter_class = argparse.RawTextHelpFormatter, 
                                        description     = "")
    parser.add_argument('--fastq1', required = True, help = "fastq File input ")
    parser.add_argument('--fastq2', required = True, help = "fastq File input ")
    parser.add_argument('--output', required = True, help = "output directory path")
    args = parser.parse_args()

    rewriteFiles(left   = args.fastq1,
                 right  = args.fastq2,
                 output = args.output)

if __name__ == "__main__":
    main()

