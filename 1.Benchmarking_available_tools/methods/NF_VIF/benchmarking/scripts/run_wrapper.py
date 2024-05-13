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

'''
This script is a wrapper for running nf-VIF
    nf-VIF only returns 3 viruses at a time in its output.
    so run nf-VIF on viruses seperately 
'''








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





    

    for virus in args.viruses.split(","):


        # create output tmp file 
        output_dir = f"{args.output}/{virus}"
        if not os.path.exists(output_dir):
            os.mkdir(output_dir)

        left = f"{output_dir}/reads_R1.fq"
        right = f"{output_dir}/reads_R2.fq"

        #~~~~~~~~~~~~~~~~~~~~~
        # Split fastq files 
        #~~~~~~~~~~~~~~~~~~~~~
        # Subset the fastq files to only include the virus of itnerest 
        logger.info(f"Virus of interest:  {virus}")
        logger.info(f"\t\t Writing left file:  {left}")
        cmd = f"./subsetHPV.py --fastq {args.fastq1} --virus {virus} --output {left}"
        subprocess.run(cmd, shell=True)

        logger.info(f"\t\t Writing right file:  {right}")
        cmd = f"./subsetHPV.py --fastq {args.fastq2} --virus {virus} --output {right}"
        subprocess.run(cmd, shell=True)






        # path to the input file
        input_file_path = f"{output_dir}/run_command.sh"
        logger.info(f"\t\t Writing command file:  {input_file_path}")

        cmd = f'''
        nextflow run /home/mbrown/GitHub/NF_VIF/nf-VIF/main.nf \
            --reads '{output_dir}/*_R{{1,2}}.fq' \
            --genome 'hg38' \
            --bwt2_index '/home/mbrown/CTAT_VIF/Other_Pipelines/VirusFiner2/Resources/human_reference' \
            --fasta   '/home/bhaas/CTAT_GENOME_LIBS/GRCh38_gencode_v22_CTAT_lib_Mar012021.plug-n-play/ctat_genome_lib_build_dir/ref_genome.fa' \
            --blatdb  '/home/mbrown/CTAT_VIF/Other_Pipelines/VirusFiner2/Resources/human_reference/GRCh38.genome.2bit' \
            --fasta_hpv '/home/mbrown/References/HPV/HPVs_db.fasta' \
            --bwt2_index_hpv '/home/mbrown/References/HPV/HPV_BWT2_index' \
            --outdir '{output_dir}' \
            -profile singularity \
            --split_report \
            --skip_trimming \
            --skip_fastqc \
            --skip_multiqc \
            --nb_geno 1
        '''

        ## Write the new input.json file 
        input_file = open(f"{output_dir}/run_command.sh", "w")
        input_file.write("#!/bin/bash\nset -e\n")
        input_file.write(cmd)
        input_file.close()


if __name__ == "__main__":
    main()











