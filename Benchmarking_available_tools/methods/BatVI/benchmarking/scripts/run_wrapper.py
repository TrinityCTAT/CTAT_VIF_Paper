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










def main():
    #add options to inputs
    parser = argparse.ArgumentParser(   formatter_class = argparse.RawTextHelpFormatter, 
                                        description     = "")
    parser.add_argument('--fastq1', required = True, help = "fastq File input ")
    parser.add_argument('--fastq2', required = True, help = "fastq File input ")
    parser.add_argument('--viruses', required = True, type=str,
                        help = "Viruses of interest")
    parser.add_argument('--output', required = True, help = "output directory path")
    args = parser.parse_args()

    # create output tmp file 
    if not os.path.exists(f"{args.output}/tmp"):
        os.mkdir(f"{args.output}/tmp")



    

    for virus in args.viruses.split(","):


        # create output tmp file 
        output_dir = f"{args.output}/{virus}"
        if not os.path.exists(output_dir):
            os.mkdir(output_dir)

        left = f"{output_dir}/left.fq"
        right = f"{output_dir}/right.fq"


        logger.info(f"Virus of interest:  {virus}")
        logger.info(f"\t\t Writing left file:  {left}")
        cmd = f"./subsetHPV.py --fastq {args.fastq1} --virus {virus} --output {left}"
        subprocess.run(cmd, shell=True)

        logger.info(f"\t\t Writing right file:  {right}")
        cmd = f"./subsetHPV.py --fastq {args.fastq2} --virus {virus} --output {right}"
        subprocess.run(cmd, shell=True)


        #~~~~~~~~~~~~~~~~~~~~~~
        # Adjust thee read names 
        #~~~~~~~~~~~~~~~~~~~~~~
        logger.info(f"\t\t Adjusting read names")
        cmd = f"./adjust_read_names.py --fastq1 {output_dir}/left.fq --fastq2 {output_dir}/right.fq --output {output_dir}"
        subprocess.run(cmd, shell=True)




        # path to the input file
        input_file_path = f"{args.output}/tmp/inputs.json"
        logger.info(f"\t\t Writing inputs file:  {input_file_path}")

        left = f"{output_dir}/noCHR_left.fq"
        right = f"{output_dir}/noCHR_right.fq"

        cmd = f'''
        "BatVI.docker":"brownmp/batvi:devel",
        "BatVI.left":"{left}",
        "BatVI.right":"{right}",
        "BatVI.Human_BLAST_index":"/home/mbrown/References/HG38/BLAST_index_noCHR.tar.gz",
        "BatVI.Human_BatIndex":"/home/mbrown/CTAT_VIF/Other_Pipelines/BatVI/Human_BatIndex.tar.gz",
        "BatVI.Human_BWA_Index":"/home/mbrown/References/HG38/BWA_index_noCHR.tar.gz",
        "BatVI.Human_fasta":"ref_genome_noCHR.fa",
        "BatVI.Virus_BLAST_index":"/home/mbrown/References/HPV/HPV_BLAST_index.tar.gz",
        "BatVI.Virus_BatIndex":"/home/mbrown/CTAT_VIF/Other_Pipelines/BatVI/Virus_BatIndex.tar.gz",
        "BatVI.Virus_BWA_Index":"/home/mbrown/References/HPV/HPV_BWA_index.tar.gz",
        "BatVI.Virus_fasta":"HPVs_db.fasta",
        "BatVI.insertion_length":"200",
        "BatVI.cpus":"1",
        "BatVI.preemptible":"1",
        "BatVI.sample_id":"TEST"
'''
        ## Add the brackets 
        cmd = "{"+cmd+"}"

        ## Write the new input.json file 
        input_file = open(f"{output_dir}/inputs.json", "w")
        input_file.write(cmd)
        input_file.close()



        # Comands to run the WDL
        cmd = """
        java -jar /home/mbrown/Tools/cromwell-71.jar \
            run /home/mbrown/CTAT_VIF/Other_Pipelines/BatVI/WDL/BatVI.wdl \
            -i inputs.json
        """

        logger.info(f"\t\t Writing command file: {output_dir}/run_command.sh")
        ## Write the new input.json file 
        ## more used for future reference 
        input_file = open(f"{output_dir}/run_command.sh", "w")
        input_file.write("#!/bin/bash\nset -e\n")
        input_file.write(cmd)
        input_file.close()

if __name__ == "__main__":
    main()











