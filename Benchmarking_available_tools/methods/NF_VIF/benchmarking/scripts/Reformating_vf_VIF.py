#!/usr/bin/env python3

import pandas as pd
import pysam 
import os, sys, re
import logging
import sys, time
import itertools
import argparse





logging.basicConfig(format='\n %(levelname)s : %(message)s', level=logging.DEBUG)
logger = logging.getLogger(__name__)

class ReformatNF_VIF:
    '''
    Class to reformat VirusFinder2 output for benchmarking 
    '''
    # initialize object
    def __init__(   self,
                    output_directory,
                    sample_tag,
                    output): # arguments to class instantiation 
        
        self.output_directory = output_directory
        self.sample_tag = sample_tag 
        self.output     = output


    def readFiles(self):
        message_str = f"\t Reading input files"
        print(message_str)

        # Get the files in the output directory 
        file_list = os.listdir(self.output_directory)
        # select the filtered outputs for each virus 
        # r = re.compile(".*_table_filtered.csv")
        r = re.compile(".*_table.csv")
        newlist = list(filter(r.match, file_list))
        df = pd.DataFrame()
        for table_file in newlist:
            message_str = f"\t\t file: {table_file}"
            print(message_str)
            output_file = os.path.join(self.output_directory, table_file )
            tmp = pd.read_csv(output_file)
            df = pd.concat([df,tmp])
        self.combined_output = df
        return self 


    def reformatOutput( self ):

        message_str = f"\t Reformating inputs"
        print(message_str)

        df = self.combined_output.rename({"genotype": "chrB", 
                                        "chr":"chrA",
                                        "chr_position":"coordA",
                                        "position":"coordB",
                                        "feature":"orientA",
                                        "strand":"orientB",
                                        "counts":"total"}, axis="columns")
        #~~~~~~~~~~
        # Adjust the orientations 
        #~~~~~~~~~~
        df['orientA'] = df['orientA'].replace(['left'],'-')
        df['orientA'] = df['orientA'].replace(['right'],'+')
        #~~~~~~~~~~
        ## Create sample name 
        #~~~~~~~~~~
        df.insert(0, "#sample", [self.sample_tag] * df.shape[0])
        # drop uneeded columns 
        # df = df.drop( labels= ["sample","multiqc_index"], axis=1)
        df = df.drop( labels= ["sample"], axis=1)
        
        # entry : create and add the entry column
        entry = df["chrA"] + "~" + df["coordA"].astype("str") + "~" + df["orientA"] + "~" + df["chrB"] + "~" + df["coordB"].astype("str") + "~" + df["orientB"]
        df.insert(1, "entry", entry)

        virus_brkend_grp = df["chrB"] + ":" + df["coordB"].astype("str")
        df.insert(10, "virus_brkend_grp", virus_brkend_grp)

        self.df = df

        return self
    
    def saveOutput( self ):
        # Save the dataframe to file 
        message_str = f"\t Saving reformatted file: {self.output}"
        print(message_str)
        self.df.to_csv( self.output, 
                       sep = "\t",
                       index = False)


def main():

    ####################
    # Parse the use supplied information 
    ####################
    # Set the variables to supply 
    parser = argparse.ArgumentParser(description="Reformat  vf-VIF output for benchamrking.", 
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--output_directory", type=str, required=False, help="output directory.", default = ".")
    parser.add_argument("--sample_tag", type=str, required=True, help="Virus identified.")
    parser.add_argument("--output", type=str, required=False, help="output file name.", default = ".")

    # Parse the variables given 
    args = parser.parse_args()
    output_directory         = args.output_directory
    output_file              = args.output
    sample_tag               = args.sample_tag

    # asign the output file if given "."
    if output_file ==  ".":
        output_file = os.getcwd()
        output_file = os.path.join(output_file,"output.tsv")


    message_str = "\n####################################################################################\n\t\t\t\tRunning Reformating\n####################################################################################"
    print(message_str)

    ##############################
    # Create object and Load Data
    ##############################
    # initiate the Vf_VIF object 
    VF = ReformatNF_VIF(output_directory = output_directory,
                        sample_tag = sample_tag,
                        output = output_file)
    VF = VF.readFiles()
    # reformat 
    VF = VF.reformatOutput()
    # Save output 
    VF.saveOutput()





if __name__ == "__main__":
    main()