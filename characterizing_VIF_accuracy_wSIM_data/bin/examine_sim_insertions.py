#!/usr/bin/env python

import sys, os, re
import argparse
import logging
import pandas
import gzip
from collections import defaultdict
from math import sqrt

logging.basicConfig(level=logging.INFO, 
                    format='%(asctime)s : %(levelname)s : %(message)s',
                    datefmt='%H:%M:%S')
logger = logging.getLogger(__name__)




def main():
    parser = argparse.ArgumentParser(description="simulation insertion site analyzer", formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("--truth_insertions", type=str, required=True, help="truth_insertions")
    parser.add_argument("--predicted_insertions", type=str, required=True, help="predicted insertions")
    parser.add_argument("--score_vir_brkend_grp_only", action='store_true', default=False, help="make truth assignments based on virus brkp end positions, but boost based on human chr distance")
    
    # other params that influence sens / specificity
    
    args = parser.parse_args()
    

    truth_filename = args.truth_insertions
    pred_insertions_filename = args.predicted_insertions
    score_vir_brkend_grp_only = args.score_vir_brkend_grp_only
    
    logger.info("-parsing truth data: {}".format(truth_filename))
    truth_data = pandas.read_csv(truth_filename, sep="\t", names=['group', 'insertion'])
    truth_data[['chrA', 'coordA_lend', 'coordA_rend', 'orientA', 'chrB', 'coordB_lend', 'coordB_rend', 'orientB']] = truth_data['insertion'].str.split("~", 7, expand=True)

    logger.info("-parsing predicted insertions: {}".format(pred_insertions_filename))
    pred_data = pandas.read_csv(pred_insertions_filename, sep="\t")

    # for the refined data
    if 'contig' in pred_data.columns:
        pred_data.rename({ 'contig' : 'entry'},
                         axis=1, inplace=True)

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
            
            insertion = Insertion(row['entry'],
                                  row['chrA'],
                                  int(row['coordA']),
                                  int(row['coordA']),
                                  row['orientA'],
                                  row['chrB'],
                                  int(row['coordB']),
                                  int(row['coordB']),
                                  row['orientB'],
                                  row)

            group_to_pred_insertions[group].append(insertion)


    # header for tab output:
    print("\t".join(["group", "truth_insertion_name",
                     "ref_chr", "truth_ref_pos", "vir_chr", "truth_vir_pos",
                     "pred_insertion_name", "pred_chr_pos", "pred_vir_pos",
                     "total_reads", "virus_brkend_grp", "is_primary"]))
    

    # perform analysis by group
    for group in group_to_truth_insertions:
        
        logger.info("-analyzing group: {}".format(group))
        
        truth_insertions = group_to_truth_insertions[group]
        pred_insertions = group_to_pred_insertions[group]
        
        analyze_insertions(group, truth_insertions, pred_insertions, score_vir_brkend_grp_only)
        


    sys.exit(0)

        
        
class Insertion:

    def __init__(self,
                 insertion_name,
                 chrA, lendA, rendA, orientA,
                 chrB, lendB, rendB, orientB,
                 row
                 ):

        self.name = insertion_name
        self.row = row
        
        # make chrA correspond to the ref genome
        if re.match("chr", chrA) and not re.match("chr", chrB):
            self.ref_chr = chrA
            self.ref_pos = rendA if orientA == '+' else lendA
            self.vir_chr = chrB
            self.vir_pos = lendB if orientB == '+' else rendB

        elif re.match("chr", chrB) and not re.match("chr", chrA):
            self.ref_chr = chrB
            self.ref_pos = lendB if orientA == '+' else rendB
            self.vir_chr = chrA
            self.vir_pos = rendA if orientB == '+' else lendA
            
        else:
            raise RuntimeError("Error, cannot resolve ref genome chr from entry: {}".format(row))



def analyze_insertions(group, truth_insertions, pred_insertions, score_vir_brkend_grp_only):

    truth_insertion_to_closest_preds = defaultdict(list)

    # map each pred insertion to the closest truth insertion and generate report
    for pred_insertion in pred_insertions:
        closest_dist = None
        closest_insertion = None
    
        for truth_insertion in truth_insertions:

            distance = None
            
            if (score_vir_brkend_grp_only == False and
                pred_insertion.ref_chr == truth_insertion.ref_chr and
                pred_insertion.vir_chr == truth_insertion.vir_chr):
                
                distance = sqrt( (pred_insertion.ref_pos - truth_insertion.ref_pos)**2 + (pred_insertion.vir_pos - truth_insertion.vir_pos)**2)


            elif (score_vir_brkend_grp_only == True and
                  pred_insertion.vir_chr == truth_insertion.vir_chr):

                distance = abs(pred_insertion.vir_pos - truth_insertion.vir_pos)
                # discount distance based on same chr
                if (pred_insertion.ref_chr == truth_insertion.ref_chr):
                    chr_distance = abs(pred_insertion.ref_pos - truth_insertion.ref_pos)
                    FUZZ_BONUS = 5
                    distance = max(0, distance - FUZZ_BONUS) - ( 1 / (chr_distance + 1) )
                    
            if distance is not None:
                if closest_insertion is None or closest_dist > distance:
                    closest_insertion = truth_insertion
                    closest_dist = distance

                    
        if closest_insertion is not None:
            truth_insertion_to_closest_preds[closest_insertion.name].append(pred_insertion)
            print("\t".join([group,
                             closest_insertion.name,
                             closest_insertion.ref_chr, str(closest_insertion.ref_pos),
                             closest_insertion.vir_chr, str(closest_insertion.vir_pos),
                             pred_insertion.name, str(pred_insertion.ref_pos), str(pred_insertion.vir_pos),
                             str(pred_insertion.row['total']) if str(pred_insertion.row['total']) != "NA" else str(pred_insertion.row['prelim.adj_total']),
                             pred_insertion.row['virus_brkend_grp'], str(pred_insertion.row['is_primary'])
                             ]
                             ))
            
        else:
            print("\t".join([group,
                             "NA",
                             pred_insertion.ref_chr, "NA",
                             pred_insertion.vir_chr, "NA",
                             pred_insertion.name, str(pred_insertion.ref_pos), str(pred_insertion.vir_pos),
                             str(pred_insertion.row['total']) if str(pred_insertion.row['total']) != "NA" else str(pred_insertion.row['prelim.adj_total']),
                             pred_insertion.row['virus_brkend_grp'], str(pred_insertion.row['is_primary'])
                             ]
                             ))
    

    # report remaining insertions that lack predictions assigned.
    for truth_insertion in truth_insertions:
        if truth_insertion.name not in truth_insertion_to_closest_preds:
            print("\t".join([group,
                             truth_insertion.name,
                             truth_insertion.ref_chr, str(truth_insertion.ref_pos),
                             truth_insertion.vir_chr, str(truth_insertion.vir_pos),
                             "NA", "NA", "NA", "NA",
                             "NA", "NA"
                             ]
                             ))
            

    return

            
            
if __name__=='__main__':
    main()



    
