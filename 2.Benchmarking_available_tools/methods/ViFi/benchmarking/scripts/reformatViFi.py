#!/usr/bin/env python3

import pandas as pd
import pysam 
import os, sys, re
import logging
import sys, time
import itertools
import numpy as np
import argparse

## Set up the logging  
logging.basicConfig(format='\n %(levelname)s : %(message)s', level=logging.DEBUG)
logger = logging.getLogger(__name__)

def table(a):
    # Function similar to R's table function 
    # Get the counts for each element present 
    dic = {}
    for i in a:
        if i in dic:
            dic[i] += 1
        else:
            dic[i] = 1
    return dic
    

class ReformatViFi:
    '''
    Class to reformat VirusFinder2 output for benchmarking 
    '''
    # initialize object
    def __init__(   self, 
                    input_file,
                    output,
                    sample_tag = "Virus_db_Dec092021.1"): # arguments to class instantiation 
        
        self.input_file = input_file
        self.sample_tag = sample_tag 
        self.output     = output

    def parseOutput( self ):

        #################
        # Read in and parse the file 
        #################
        message_str = f"\t\tReading input file {self.input_file}"
        logger.info(message_str)

        sep_string = "##==========================================================================================================================================================================================================================\n"
        # Seperate all individual insertions into lists 
        insertion_list = []
        new_list = []
        counter=0
        with open(self.input_file) as fileIn:
            for i in fileIn:
                if i == sep_string:
                    counter += 1
                    insertion_list.append(new_list)
                    new_list = []
                else:
                    new_list.append(i)
            # need to add the last one 
            insertion_list.append(new_list)
        # delete the first entry which is the header 
        del(insertion_list[0])

        ######################
        # Get the first line of each insertions. The insertion info 
        ######################
        insertion_info = []
        for i in insertion_list:
            insertion_info.append(i.pop(0).rstrip().split("\t"))
        #make sure have the correct amoutn 
        #########
        #CHECK
        print(counter == len(insertion_info))

        self.insertion_info  =  insertion_info




        #############################
        # Check what reads are mapped to 
        #############################
        # ALL of the insertions
        all_insertions = []
        counter = 0
        for i in insertion_list:
            all_insertions.append([])
            for j in i:
                a = j.split("\t")[1]
                all_insertions[counter].append(a)
            counter+=1
        # all
        
        # # table counts for each insertion 
        # for i in all_insertions:
        #     print(table(i))

        self.all_insertions = all_insertions

        return self 


    def reformat( self ):

        message_str = f"\t\tReformating the file."
        logger.info(message_str)

        #####################
        # Put into dataframe
        #####################
        df = pd.DataFrame(self.insertion_info, columns= ["chrA","Min","Max","total_reads", "L_reads","R_reads", "min","max","split1","split2"]) 
        ids = df["chrA"] + ":" + df["Min"] + "-" + df["Max"]

        # create entry field
        ids = df[["chrA", "Min","Max"]].agg('~'.join, axis=1)
        df.insert(0, "entry", ids)

        #~~~~~~~~~~~~~~~~~~~~~
        # Create CordA
        #~~~~~~~~~~~~~~~~~~~~~
        # Contert to int
        df[["Min","Max","total_reads", "L_reads","R_reads", "min","max","split1","split2"]] = df[["Min","Max","total_reads", "L_reads","R_reads", "min","max","split1","split2"]].astype(int)
        df = df.replace({"split1" : {-1 : np.nan}, 
                          "split2" : {-1 : np.nan}})
        # Replace na with 0's 
        df = df.fillna(0)

        coordA = df["split1"] + df["split2"]
        df.insert(2, "coordA", coordA)
        df

        #~~~~~~~~~~~~~~~~~~~~~
        # Generate Orientation
        #~~~~~~~~~~~~~~~~~~~~~
        orientation = []
        for i in list(zip(df["L_reads"],df["R_reads"])):
            if i[0] == 0:
                orientation.append("-")
            else: 
                orientation.append("+")
        # insert into df 
        df.insert(3, "orientA", orientation)




        ######################################
        # Get the virus insertion type 
        ######################################
        # table counts for each insertion 
        # Get the top 2 highest values, these are our contigs 
        def sortTuple(tup_list):
            # Sort values 
            ## sort tuple by second value 
            tup_list.sort(key = lambda x: x[1], reverse=True)
            return tup_list

        values = []
        for i in self.all_insertions:
            # Sort values 
            tup = list(zip(table(i), table(i).values()))
            ## get the top two 
            top_two = sortTuple(tup)[0:2]
            values.append(top_two)
        virus_type = [i[1][0] for i in values]
        virus_total = [i[1][1] for i in values]

        df.insert(4, "chrB", virus_type)
        df.insert(5, "total", virus_total)


        #~~~~~~~~~~~~~~~~~~~~
        # Change column names 
        #~~~~~~~~~~~~~~~~~~~~
        df = df.rename(columns = {"Min" : "coordA_lend", "Max" : "coordA_rend"})

        #~~~~~~~~~~~~~~~~~~~~
        # Add Sample Tag
        #~~~~~~~~~~~~~~~~~~~~
        samples = [self.sample_tag] * len(virus_total)
        df.insert(0, "#sample", samples)
        # final_df = df[df.chrB.isin(["HPV45","HPV16","HPV18"])]

        self.df  = df
        return self 

    def saveOutput( self ):
        output_file = "ViFi_Output_reformated.txt"
        output_path_str = os.path.join(self.output, output_file)
        logger.info(f"\t\t Outputing file to: {output_path_str}")

        self.df.to_csv( output_path_str, 
                        sep = "\t", index = False)
        





def main():

    ####################
    # Parse the use supplied information 
    ####################
    # Set the variables to supply 
    parser = argparse.ArgumentParser(description="Reformat the output given by ViFi to  allow it to be benchmarked.", 
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--inputs", type=str, required=True, help="truth_insertions")
    parser.add_argument("--sampleID", type=str, required=True, help="Sample id that resembles the simulated group. ex) Virus_db_Dec092021.1")
    # parser.add_argument("--Viruses", type=str, required=False, help="Viruses to consider seperated by commas.", default = "HPV45,HPV16,HPV18")
    parser.add_argument("--output", type=str, required=False, help="output directory.", default = ".")

    # Parse the variables given 
    args = parser.parse_args()
    inputs = args.inputs
    sampleID = args.sampleID
    # Viruses = args.Viruses.split(",")
    output = args.output

    if output ==  ".":
        output = os.getcwd()

    # initiate the ViFi object 
    vifi = ReformatViFi( input_file = inputs,
                         output     = output,
                         sample_tag = sampleID)
    vifi = vifi.parseOutput()
    vifi =  vifi.reformat()
    vifi = vifi.saveOutput()











if __name__ == "__main__":
    main()



















