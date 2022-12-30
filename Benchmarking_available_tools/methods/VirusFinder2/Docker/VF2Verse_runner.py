#!/usr/bin/env python

import sys, os, re
import subprocess
import argparse
import logging


logging.basicConfig(stream=sys.stderr, level=logging.INFO)
logger = logging.getLogger(__name__)

utildir = os.path.abspath(os.path.dirname(__file__))


def main():

    parser = argparse.ArgumentParser(description="run VF2/Verse", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--human_resource_dir", type=str, required=True, help="human resource dir path")
    parser.add_argument("--virus_resource_dir", type=str, required=True, help="virus resource dir path")
    parser.add_argument("--target_virus", type=str, required=True, help="virus to examine")
    parser.add_argument("--left_fq", type=str, required=True, help="left fastq")
    parser.add_argument("--right_fq", type=str, required=True, help="right fastq")
    parser.add_argument("--sample_id", type=str, required=True, help="sample id")
    parser.add_argument("--CPU", type=int, default=5, help='cpu count for multithreading')

    
    args = parser.parse_args()


    cmd = f"{utildir}/write_configuration_file.py  --fastq1 {args.left_fq} --fastq2 {args.right_fq} --human_resource_dir {args.human_resource_dir} --virus_resource_dir {args.virus_resource_dir}"
    logger.info(cmd)
    subprocess.check_call(cmd, shell=True)


    cmd = f"{utildir}/VirusFinder2.0/preprocess.pl -c configuration.txt"
    logger.info(cmd)
    subprocess.check_call(cmd, shell=True)


    cmd = f"{utildir}/VirusFinder2.0/detect_integration.pl -c configuration.txt -v {args.target_virus}"
    logger.info(cmd)
    subprocess.check_call(cmd, shell=True)

    assert os.path.exists("results-virus-loci.txt"), "Error - missing results output file 'results-virus-loci.txt'"
    
    os.rename("results-virus-loci.txt", f"{args.target_virus}.VF2Verse.txt")
    
    logger.info(f"Done, see output file: {args.target_virus}.VF2Verse.txt")
    
    sys.exit(0)
    


if __name__=='__main__':
    main()
