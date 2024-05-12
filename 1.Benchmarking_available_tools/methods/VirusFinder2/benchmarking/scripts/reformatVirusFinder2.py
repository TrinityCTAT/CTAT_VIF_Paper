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

class ReformatVirusFinder2:
    '''
    Class to reformat VirusFinder2 output for benchmarking 
    '''
    # initialize object
    def __init__(   self,
                    input_file,
                    virus_type,
                    sample_tag,
                    output): # arguments to class instantiation 
        
        self.input_file = input_file
        self.virus_type = virus_type
        self.sample_tag = sample_tag 
        self.output_path = output 

    def reformatOutput( self ):

        ## read in the input file 
        df = pd.read_csv(self.input_file,
                         sep="\t")
        df["Chromosome 1"] = df["Chromosome 1"].replace({"chrVirus":self.virus_type})
        df["Chromosome 2"] = df["Chromosome 2"].replace({"chrVirus":self.virus_type})
        header = ["chrA","coordA","orientA", "chrB", "coordB", "orientB", "#Support reads (pair+softclip)", "Confidence"]
        df.columns = header
        df

        ## Create sample name 
        df.insert(0, "#sample", [self.sample_tag] * df.shape[0])

        ## Create contig Name
        ### adjust orientB if not given 
        idx = df["orientB"].isnull()
        df.loc[idx,["orientB"]] = "NA"
        ### combine 
        contig = df["chrA"] + "~" + df["coordA"].astype(str) + "~" + df["orientA"] + "~" + df["chrB"] + "~" + df["coordB"].astype(str) + "~" + df["orientB"]
        df.insert(1,"entry",contig)

        ## get the total number of suporting  reads 
        sums = df['#Support reads (pair+softclip)'].str.split("+", expand=True).astype(int).sum(axis =1)
        df["total"] = sums


        idx = df.chrB.str.match("^chr*")
        df.loc[idx, ['chrA','coordA','orientA','chrB','coordB','orientB']] = df.loc[idx, ['chrB','coordB','orientB','chrA','coordA','orientA']].values
        
        # Print nice output for refereencee
        tmp = df.iloc[:5].to_string()
        tmp = tmp.split("\n")
        tmp = [f"\t\t\t\t{i}\n" for i in tmp]
        tmp = "".join(tmp)
        print(f"\t Examine output:\n{tmp}")

        self.df = df

        return self

    def saveOutput( self ):

        ## Save the output file as a tsv 
        output_file = f"{self.virus_type}_VirusFindeer2_output.tsv" 
        logger.info(f"\t\tOutput File: {output_file}")

        output_file = os.path.join(self.output_path, output_file)

        ## SAVE
        self.df.to_csv(output_file, sep = "\t", index = False)





def main():

    ####################
    # Parse the use supplied information 
    ####################
    # Set the variables to supply 
    parser = argparse.ArgumentParser(description="Reformat  VirusFinder2 output for benchamrking.", 
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--predicted_insertions", type=str, required=True, help="predicted insertions")
    parser.add_argument("--virus", type=str, required=True, help="Virus identified.")
    parser.add_argument("--sample_tag", type=str, required=True, help="Virus identified.")
    parser.add_argument("--output", type=str, required=False, help="output directory.", default = ".")

    # Parse the variables given 
    args = parser.parse_args()
    pred_insertions_filename = args.predicted_insertions
    virus = args.virus
    output_path = args.output
    sample_tag = args.sample_tag

    if output_path ==  ".":
        output_path = os.getcwd()

    message_str = "\n####################################################################################\n\t\t\t\tRunning Reformating\n####################################################################################"
    print(message_str)

    ##############################
    # Load Data
    ##############################
    # initiate the ViFi object 
    VirusFinder = ReformatVirusFinder2( input_file = pred_insertions_filename,
                                        virus_type = virus,
                                        sample_tag = sample_tag,
                                        output     = output_path)
    
    VirusFinder = VirusFinder.reformatOutput()

    VirusFinder.saveOutput()




if __name__ == "__main__":
    main()



