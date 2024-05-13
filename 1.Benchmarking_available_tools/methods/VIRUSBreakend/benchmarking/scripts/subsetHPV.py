#!/usr/bin/env python3




import pandas as pd
import pysam 
import os, sys, re
import logging
import sys, time
import itertools
import argparse
logging.basicConfig(format='\n %(levelname)s : %(message)s', level=logging.DEBUG)


'''
The simulated data has virus insertions. 
Use this script to pull out the alignments with the simulated virus of interest 
'''


exec(open('FastFileClass.py').read())




def subsetHPV(sequences_dic, virus):
    a = [i for i in sequences_dic]
    b = [i.split("~")[4] for i in a]
    c = [i == virus for i in b]
    d = [i for i, x in enumerate(c) if x]
    return [a[i] for i in d]


def main():
    #add options to inputs
    parser = argparse.ArgumentParser(   formatter_class = argparse.RawTextHelpFormatter, 
                                        description     = "")
    parser.add_argument('--fastq', required = True, help = "fastq File input ")
    parser.add_argument('--virus', required = True, help = "Virus of interest")
    parser.add_argument('--output', required = True, help = "output directory and file path")
    args = parser.parse_args()

    inputfile = fastFile(input_file = args.fastq)

    # Read in the file 
    inputfile = inputfile.fqReader()

    new_ids = subsetHPV(  sequences_dic = inputfile.sequences, 
                                virus = args.virus)

    new_sequences = {i:inputfile.sequences[i] for i in new_ids}

    inputfile.sequences = new_sequences

    inputfile.fqWriter(args.output)

if __name__ == "__main__":
    main()















