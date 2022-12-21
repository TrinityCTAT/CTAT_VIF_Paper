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
logging.basicConfig(format="\n %(levelname)s : %(message)s", level=logging.DEBUG)
logger = logging.getLogger(__name__)


class ReformatBatVI:
    """
    Class to reformat BatVI output for benchmarking 
    """

    # initialize object
    def __init__(
        self, input_file, Virus, output, sample_tag="test"
    ):  # arguments to class instantiation

        self.input_file = input_file
        self.Virus = Virus
        self.sample_tag = sample_tag
        self.output = output

    def reformat(self):

        message_str = f"\t\tReformating the file."
        logger.info(message_str)

        #####################
        # Put into dataframe
        #####################
        # column_names = ["chrA","Min","Max","total_reads", "L_reads","R_reads", "min","max","split1","split2"]
        df = pd.read_csv(self.input_file, sep="\t", index_col=False)

        samples = [self.sample_tag] * df.shape[0]
        df.insert(1, "#sample", samples)
        df = df.drop(labels=["LIB"], axis=1)

        # Rename the columns
        df = df.rename(
            {
                "Chr": "chrA",
                "Human Pos": "coordA",
                "Sign": "orientA",
                "Viral Sign": "orientB",
                "Viral Pos": "coordB",
                "Read Count": "total",

                # replace spaceswith underscores
                "Split Reads" : "Split_Reads",
                "Uniquely Mapped Reads", "Uniquely_Mapped_Reads",
                "Multiply Mapped Reads", "Multiply_Mapped_Reads",
                "Rank1 Hits", "Rank1_Hits"
            },
            axis="columns",
        )

        # chrA : adjust the format, add 'chr' to begining
        df["chrA"] = "chr" + df["chrA"].astype(str)

        # Add virus of interest
        virus_list = [self.Virus] * df.shape[0]
        df.insert(5, "chrB", virus_list)

        # entry : create and add the entry column
        entry = (
            df["chrA"]
            + "~"
            + df["coordA"].astype("str")
            + "~"
            + df["chrB"]
            + "~"
            + df["coordB"].astype("str")
        )
        df.insert(1, "entry", entry)

        # Reorder columns
        column_order = [
            "#sample",
            "entry",
            "chrA",
            "coordA",
            "orientA",
            "chrB",
            "coordB",
            "orientB",
            "total",
            "Split_Reads",
            "Uniquely_Mapped_Reads",
            "Multiply_Mapped_Reads",
            "Rank1_Hits",
        ]
        df = df[column_order]

        self.df = df
        return self

    def saveOutput(self):

        ## Save the output file as a tsv
        output_file = f"{self.Virus}_BatVI_output.tsv"
        logger.info(f"\t\tOutput File: {output_file}")

        ## SAVE
        self.df.to_csv(output_file, sep="\t", index=False)


def main():

    ####################
    # Parse the use supplied information
    ####################
    # Set the variables to supply
    parser = argparse.ArgumentParser(
        description="Reformat  BatVI output for benchamrking.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--predicted_insertions", type=str, required=True, help="predicted insertions"
    )
    parser.add_argument("--virus", type=str, required=True, help="Virus identified.")
    parser.add_argument(
        "--sample_tag", type=str, required=True, help="Virus identified."
    )
    parser.add_argument(
        "--output", type=str, required=False, help="output directory.", default="."
    )

    # Parse the variables given
    args = parser.parse_args()
    pred_insertions_filename = args.predicted_insertions
    virus = args.virus
    output_path = args.output
    sample_tag = args.sample_tag

    if output_path == ".":
        output_path = os.getcwd()

    message_str = "\n####################################################################################\n\t\t\t\tRunning Reformating\n####################################################################################"
    print(message_str)

    ##############################
    # Load Data
    ##############################
    # initiate the ViFi object
    BatVI = ReformatBatVI(
        input_file=pred_insertions_filename,
        Virus=virus,
        sample_tag=sample_tag,
        output=output_path,
    )

    BatVI = BatVI.reformat()

    BatVI.saveOutput()


if __name__ == "__main__":
    main()
