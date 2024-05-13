#!/usr/bin/env python3

import pandas as pd
import pysam 
import os, sys, re
import logging
import sys, time
import itertools
import argparse
import subprocess


logging.basicConfig(format='\n %(levelname)s : %(message)s', level=logging.DEBUG)
logger = logging.getLogger(__name__)

utildir = os.path.abspath(os.path.dirname(__file__))


class benchmarkWrapper:
    '''
    Class to reformat VirusFinder2 output for benchmarking 
    '''
    # initialize object
    def __init__(   self,
                    output_directory,
                    viruses,
                    output): # arguments to class instantiation 
        
        self.output_directory = output_directory
        self.viruses          = viruses.split(",")
        self.output           = output
        self.out_dirs         = []


    def fileList(self):
        message_str = f"\tgetting directories"
        print(message_str)
        
        # Get the files in the output directory 
        dir_list = [os.path.join(self.output_directory,i) for i in self.viruses]
        file_list = [os.path.join(self.output_directory,i,"hpv_mapping/blat") for i in self.viruses]

        out_dirs = zip(dir_list,file_list,self.viruses)

        self.out_dirs = out_dirs
        return self 


    def runReformating( self ):

        for i in self.out_dirs:
            print(f"Virus: {i[2]}")
            #cmd = f"""./Reformating_vf_VIF.py \n\t --output_directory {i[1]} \n\t --sample_tag subset.insertion_seqs.fa \n\t --output {i[0]}"""

            cmd = f"""{utildir}/Reformating_vf_VIF.py \
                        --output_directory {i[1]} \
                        --sample_tag subset.insertion_seqs.fa \
                        --output {i[0]}/vf_VIF_output.tsv"""
            print(cmd)
            subprocess.check_call(cmd, shell=True)


    def joinReformatedFiles( self ):
        '''
        Join the reformated files together 
        '''
        message_str = "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\t\t\t\tJoining Files\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
        print(message_str)

        # Concatinate the dataframes 
        df = pd.DataFrame()
        for i in self.out_dirs:
            table_file = f"{i[0]}/vf_VIF_output.tsv"
            message_str = f"\t\t file: {table_file}"
            print(message_str)

            # load file in 
            tmp = pd.read_csv(table_file, sep = "\t")
            df = pd.concat([df,tmp])
        self.df = df

        print(df)
        
        return self 


    def saveOutput( self ):
        # Save the dataframe to file 
        message_str = f"\nSaving Final Combined file: {self.output}"
        print(message_str)
        self.df.to_csv( self.output,
                        sep = "\t",
                        index = False )



def main():

    ####################
    # Parse the use supplied information 
    ####################
    # Set the variables to supply 
    parser = argparse.ArgumentParser(description="Reformat  vf-VIF output for benchamrking.", 
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--output_directory", type=str, required=False, help="output directory.", default = ".")
    parser.add_argument('--viruses', required = True, type=str,
                        help = "Viruses of interest")
    parser.add_argument("--output", type=str, required=False, help="output file name.", default = ".")

    # Parse the variables given 
    args = parser.parse_args()
    output_directory         = args.output_directory
    output_file              = args.output
    viruses                  = args.viruses

    # asign the output file if given "."
    if output_file ==  ".":
        output_file = os.getcwd()
        output_file = os.path.join(output_file,"output.tsv")


    message_str = "\n####################################################################################\n\t\t\t\tRunning Wrapper\n####################################################################################"
    print(message_str)

    ##############################
    # Create object and Load Data
    ##############################
    # initiate the Vf_VIF object 
    VF = benchmarkWrapper(output_directory = output_directory,
                        viruses = viruses,
                        output = output_file)
    VF = VF.fileList()
    # reformat 
    VF.runReformating()

    # join together
    VF = VF.fileList()
    VF = VF.joinReformatedFiles()
    VF.saveOutput()






if __name__ == "__main__":
    main()
