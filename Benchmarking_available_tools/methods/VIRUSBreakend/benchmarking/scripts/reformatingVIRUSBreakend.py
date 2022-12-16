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

class ReformatVIRUSBreakend:
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
        self.output     = output


    def readFile(self):
        file_open = open(self.input_file,  "r")
        file_input = file_open.readlines()
        file_open.close()

        new = []
        for i in file_input:
            if i[0:2] != "##":
                new.append(i.rstrip().split("\t"))
        # get the header
        header = new.pop(0)
        # create new dataframe
        df = pd.DataFrame(new)
        df.columns = header

        self.integrations = df

        return self

    def reformatOutput( self ):


        # Fix the ID column for merging
        self.integrations["ID"] = self.integrations["ID"].str.replace("_host","")
        self.integrations["ID"] = self.integrations["ID"].str.replace("_virus","")
        
        # seperate the human and virus integrations 
        # chr_idx = self.integrations["#CHROM"].str.contains("chr")
        chr_idx = self.integrations["REF"].str.contains("N")
        human_idx = [i for i,j in enumerate(chr_idx) if j == True]
        viral_idx = [i for i,j in enumerate(chr_idx) if j == False]

        human_df = self.integrations.iloc[human_idx]
        
        # combine the viral POS to the correct human integration point 
        viral_df = self.integrations.iloc[viral_idx]
        viral_df = viral_df.rename({"POS":"coordB"}, axis='columns')
        viral_df_pos = viral_df[["coordB","ID"]]
        df = human_df.merge(viral_df_pos, on="ID")

        coordB = df.pop("coordB")
        df.insert(3, "coordB", coordB)

        # Adjust column names
        df = df.rename({"#CHROM":"chrA","POS":"coordA"},
                       axis='columns')
        
        # human orientation
        orientation = [i.split("|")[1] for i in df["INFO"]]
        df.insert(3, "orientA", orientation)
        # insert virus type 
        df.insert(4, "chrB", self.virus_type)
        
        


        #~~~~~~~~~~~~~~
        # Get total reads
        #~~~~~~~~~~~~~~
        BVF_list = []
        r = re.compile("BVF=*")
        for i in viral_df["INFO"]:
            tmp = i.split(";")
            BVF = list(filter(r.match, tmp))[0]
            BVF = BVF.replace("BVF=","")
            BVF_list.append(BVF)

        df.insert(6, "total", BVF_list)
        
        #~~~~~~~~~~
        # Create entry
        #~~~~~~~~~~
        entry = df["chrA"] + "~" + df["coordA"].astype(str) + "~" + df["orientA"] + "~" + df["chrB"] + "~" + df["coordB"].astype(str)
        df.insert(0, "entry", entry)

        #~~~~~~~~~~
        ## Create sample name 
        #~~~~~~~~~~
        df.insert(0, "#sample", [self.sample_tag] * df.shape[0])


        df = df.drop( labels= ["ID"], axis=1)


        

        self.df = df

        return self


    def saveOutput( self ):

        ## Save the output file as a tsv 
        output_file = f"{self.virus_type}_VIRUSBreakend_output.tsv" 
        full_output_file = os.path.join(self.output,output_file)
        logger.info(f"\t\tOutput File: {full_output_file}")

        ## SAVE
        self.df.to_csv(full_output_file, sep = "\t", index = False)





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

    ##############################
    # Load Data
    ##############################
    # initiate the ViFi object 
    VB = ReformatVIRUSBreakend( input_file = pred_insertions_filename,
                                        virus_type = virus,
                                        sample_tag = sample_tag,
                                        output     = output_path)
    



    message_str = "\n####################################################################################\n\t\t\t\tRunning Reformating\n####################################################################################"
    print(message_str)
    VB = VB.readFile()
    VB = VB.reformatOutput()




    message_str = "\n####################################################################################\n\t\t\t\tSaving Output\n####################################################################################"
    print(message_str)
    VB.saveOutput()




if __name__ == "__main__":
    main()



