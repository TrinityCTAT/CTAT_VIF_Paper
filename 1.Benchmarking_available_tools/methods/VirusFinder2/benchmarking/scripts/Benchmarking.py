#!/usr/bin/env python3


import sys, os, re
import argparse
import logging
import pandas
import gzip
from collections import defaultdict
from math import sqrt

## Set up the logging  
logging.basicConfig(format='\n %(levelname)s : %(message)s', level=logging.DEBUG)
logger = logging.getLogger(__name__)




        
class Insertion:
    '''
    Class for insertions 
    '''
    def __init__(self,
                 insertion_name,
                 chrA, lendA, rendA, orientA,
                 chrB, lendB, rendB, orientB,
                 row
                 ):

        self.name = insertion_name
        self.row = row
        
        # make chrA correspond to the ref genome
        ## CHR is in chrA
        if re.match("chr", chrA) and not re.match("chr", chrB):
            self.ref_chr = chrA
            self.ref_pos = rendA if orientA == '+' else lendA
            self.vir_chr = chrB
            self.vir_pos = lendB if orientB == '+' else rendB
        ## CHR is in chrB
        elif re.match("chr", chrB) and not re.match("chr", chrA):
            self.ref_chr = chrB
            self.ref_pos = lendB if orientA == '+' else rendB
            self.vir_chr = chrA
            self.vir_pos = rendA if orientB == '+' else lendA
            
        else:
            raise RuntimeError("Error, cannot resolve ref genome chr from entry: {}".format(row))

class predInsertion:
    '''
    Seperate class for insertions 
    '''
    def __init__(self,
                 insertion_name,
                 chrA, lendA, rendA, orientA,
                 chrB,
                 row
                 ):

        self.name = insertion_name
        self.row = row
        
        # make chrA correspond to the ref genome
        ## CHR is in chrA
        if re.match("chr", chrA) and not re.match("chr", chrB):
            self.ref_chr = chrA
            self.ref_pos = rendA if orientA == '+' else lendA
            self.vir_chr = chrB
            
        else:
            raise RuntimeError("Error, cannot resolve ref genome chr from entry: {}".format(row))

def analyze_insertions(group, truth_insertions, pred_insertions):

    truth_insertion_to_closest_preds = defaultdict(list)

    final_output = []
    # map each pred insertion to the closest truth insertion and generate report
    for pred_insertion in pred_insertions:
        closest_dist = None
        closest_insertion = None
    
        # Compare the predicted and truth insertions 
        ## find the ones that match togetheer 
        for truth_insertion in truth_insertions:
            if ( (pred_insertion.ref_chr == truth_insertion.ref_chr) and 
                 (pred_insertion.vir_chr == truth_insertion.vir_chr) ):
                # distance = sqrt( (pred_insertion.ref_pos - truth_insertion.ref_pos)**2 + (pred_insertion.vir_pos - truth_insertion.vir_pos)**2)
                distance = sqrt( (pred_insertion.ref_pos - truth_insertion.ref_pos )**2)

                if closest_insertion is None or closest_dist > distance:
                    closest_insertion = truth_insertion
                    closest_dist = distance
                
        if closest_insertion is not None:
            truth_insertion_to_closest_preds[closest_insertion.name].append(pred_insertion)
            output = ("\t".join([group,
                             closest_insertion.name,
                             closest_insertion.ref_chr, str(closest_insertion.ref_pos),
                             closest_insertion.vir_chr, 
                             #str(closest_insertion.vir_pos),
                             pred_insertion.name, 
                             str(pred_insertion.ref_pos), 
                             #str(pred_insertion.vir_pos), 
                             str(pred_insertion.row['total']) ]))
            
        else:
            output = ("\t".join([group,
                             "NA",
                             pred_insertion.ref_chr, "NA",
                             pred_insertion.vir_chr,
                             pred_insertion.name, 
                             str(pred_insertion.ref_pos), 
                             #str(pred_insertion.vir_pos), 
                             str(pred_insertion.row['total']) ]))
        final_output.append(output)

    # report remaining insertions that lack predictions assigned.
    for truth_insertion in truth_insertions:
        if truth_insertion.name not in truth_insertion_to_closest_preds:
            output = ("\t".join([group,
                             truth_insertion.name,
                             truth_insertion.ref_chr, str(truth_insertion.ref_pos),
                             truth_insertion.vir_chr, 
                             #str(truth_insertion.vir_pos),
                            #  "NA", 
                             "NA", "NA", "NA"]))
            final_output.append(output)
            

    return final_output
    





def main():

    ####################
    # Parse the use supplied information 
    ####################
    # Set the variables to supply 
    parser = argparse.ArgumentParser(description="simulation insertion site analyzer", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--truth_insertions", type=str, required=True, help="truth_insertions")
    parser.add_argument("--predicted_insertions", type=str, required=True, help="predicted insertions")
    parser.add_argument("--Viruses", type=str, required=False, help="Viruses to consider seperated by commas.", default = None)
    parser.add_argument("--output", type=str, required=False, help="output directory.", default = ".")

    # Parse the variables given 
    args = parser.parse_args()
    truth_filename = args.truth_insertions
    pred_insertions_filename = args.predicted_insertions
    Viruses = args.Viruses
    output_path = args.output

    if output_path ==  ".":
        output_path = os.getcwd()


    message_str = "\n####################################################################################\n\t\t\t\tRunning Benchmarking\n####################################################################################"
    print(message_str)

    ##############################
    # Load Data
    ##############################
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # load in the truth set data 
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    logger.info("-parsing truth data: {}".format(truth_filename))
    truth_data = pandas.read_csv(truth_filename, sep="\t", names=['group', 'insertion'])
    truth_data[['chrA', 'coordA_lend', 'coordA_rend', 'orientA', 'chrB', 'coordB_lend', 'coordB_rend', 'orientB']] = truth_data['insertion'].str.split("~", 7, expand=True)


    #~~~~~~~~~~~~~~~~~~~~
    # Subset the truthset to  only viruses of interst 
    #~~~~~~~~~~~~~~~~~~~~
    before_shape = truth_data.shape

    if Viruses != None:
        # parse the virusees 
        Viruses = Viruses.split(",")

        print(f"\tViruses  of interest: {Viruses}")
        # Subset the truth set to only included the wanted viruses 
        truth_data = truth_data[truth_data.chrB.isin(Viruses)]

        after_shape = truth_data.shape
        print(f"\tBefore Subset:\n\t\t{before_shape[0]}\n\tAfter Subset:\n\t\t{after_shape[0]}")

    print(f"\tTruthset dimensions:\n\t\t{before_shape[0]}")
    ### Print it nicely 
    df_string = truth_data.head().to_string()
    df_string = df_string.replace("\n","\n\t")
    print("Truthset dataframe head:\n",df_string)

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Load in predicted data 
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    logger.info("-parsing predicted insertions: {}".format(pred_insertions_filename))
    pred_data = pandas.read_csv(pred_insertions_filename, sep="\t")

    
    ### Print it nicely 
    df_string = pred_data.head().to_string()
    df_string = df_string.replace("\n","\n\t")
    print("Preicted dataframe head:\n",df_string)


    group_to_truth_insertions = defaultdict(list)

    # organize by group:

    ## organize truth insertions
    for _, row in truth_data.iterrows():

        group = row['group']
        
        insertion = Insertion(row['insertion'],
                                row['chrA'],
                                int(row['coordA_lend']),
                                int(row['coordA_rend']),
                                row['orientA'],
                                row['chrB'],
                                int(row['coordB_lend']),
                                int(row['coordB_rend']),
                                row['orientB'],
                                row)

        group_to_truth_insertions[group].append(insertion)


    ## organize predicted insertions
    group_to_pred_insertions = defaultdict(list)

    for _, row in pred_data.iterrows():
        group = row['#sample']

        if row['total'] > 0 and (  (re.match("chr", row['chrA']) is not None) ^ (re.match("chr", row['chrB']) is not None) ):
            
            insertion = predInsertion(row['entry'],
                                    row['chrA'],
                                    int(row['coordA']),
                                    int(row['coordA']),
                                    row['orientA'],
                                    row['chrB'],
                                    row)

            group_to_pred_insertions[group].append(insertion)


    logger.info("Preparing output file.")
    # header for tab output:
    header = "\t".join(["group", "truth_insertion_name",
                        "ref_chr", "truth_ref_pos", "vir_chr", 
                        #"truth_vir_pos",
                        "pred_insertion_name", 
                         "pred_chr_pos", 
                        # "pred_vir_pos",
                        "total_reads"])

    final_list = [header]
    
    # perform analysis by group
    for group in group_to_truth_insertions:
        
        logger.info("-analyzing group: {}".format(group))
        
        truth_insertions = group_to_truth_insertions[group]
        pred_insertions = group_to_pred_insertions[group]
        
        output = analyze_insertions(group, truth_insertions, pred_insertions)
        final_list = final_list + output
    final_string = "\n".join(final_list)


    #########################
    # output results to the file 
    #########################
    output_file_str = os.path.join(output_path, "Benchmarking_Output.txt")
    logger.info( f"Outputing results to: {output_file_str}" )
    output_file = open(output_file_str, "w")
    output_file.write(final_string)
    output_file.close()

        

if __name__ == "__main__":
    main()

















