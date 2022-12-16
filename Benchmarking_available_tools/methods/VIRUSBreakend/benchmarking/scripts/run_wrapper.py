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


    CWD = os.getcwd()
    

    for virus in args.viruses.split(","):
        print(f"################################\n\tRunning VIRUSBreakend:  {virus}\n################################\n")


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

        #~~~~~~~~~~~~~~~~~~~~~~~
        #Run VIRUSBreakend
        #~~~~~~~~~~~~~~~~~~~~~~~
        # Chanage directory
        os.chdir(output_dir)


        logger.info(f"\t\t Running VIRUSBreakend")
        cmd = f"""java -jar /home/mbrown/Tools/cromwell-71.jar run ../WDL/VIRUSBreakend.wdl -i ../WDL/inputs.json"""

        subprocess.run(cmd, shell=True)

        # remove the memmory exhausting inputs 
        cmd ="rm -r cromwell-executions/VIRUSBreakend/*/call-RunVIRUSBreakend/inputs"
        subprocess.run(cmd, shell=True)

        # Chanage back to original directory
        os.chdir(args.output)


if __name__ == "__main__":
    main()
