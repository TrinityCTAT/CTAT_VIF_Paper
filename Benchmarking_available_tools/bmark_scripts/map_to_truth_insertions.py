#!/usr/bin/env python

import sys, os, re
import argparse
import logging
import pandas
import gzip
from collections import defaultdict
from math import sqrt
import csv

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s : %(levelname)s : %(message)s",
    datefmt="%H:%M:%S",
)
logger = logging.getLogger(__name__)


def main():
    parser = argparse.ArgumentParser(
        description="simulation insertion site analyzer",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "--truth_insertions", type=str, required=True, help="truth_insertions"
    )
    parser.add_argument(
        "--predicted_insertions", type=str, required=True, help="predicted insertions"
    )

    # other params that influence sens / specificity

    args = parser.parse_args()

    truth_filename = args.truth_insertions
    pred_insertions_filename = args.predicted_insertions

    logger.info("-parsing truth data: {}".format(truth_filename))
    truth_data = pandas.read_csv(
        truth_filename, sep="\t", names=["sample_name", "insertion"]
    )
    truth_data[
        [
            "chrA",
            "coordA_lend",
            "coordA_rend",
            "orientA",
            "chrB",
            "coordB_lend",
            "coordB_rend",
            "orientB",
        ]
    ] = truth_data["insertion"].str.split("~", 7, expand=True)

    group_to_truth_insertions = defaultdict(list)

    # organize by group:

    ## organize truth insertions
    for _, row in truth_data.iterrows():

        sample_name = row["sample_name"]
        human_chr, human_brkpt, virus, virus_brkpt = extract_truth_breakpoint_info(row)
        insertion = Insertion(human_chr, human_brkpt, virus, virus_brkpt, row)
        group_to_truth_insertions[sample_name].append(insertion)

    ################################
    ## organize predicted insertions
    group_to_pred_insertions = defaultdict(list)

    logger.info("-parsing predicted insertions: {}".format(pred_insertions_filename))
    pred_data = pandas.read_csv(pred_insertions_filename, sep="\t")

    column_names = pred_data.columns.values.tolist()

    for _, row in pred_data.iterrows():
        group = row["sample_name"]

        insertion = Insertion(
            row["human_chr"], row["human_brkpt"], row["virus"], row["virus_brkpt"], row
        )

        group_to_pred_insertions[group].append(insertion)

    column_names.insert(0, "pred_insertion_name")
    column_names.extend(
        [
            "truth_insertion_name",
            "truth_human_chr",
            "truth_human_brkpt",
            "truth_virus",
            "truth_virus_brkpt",
        ]
    )

    NA_row = {i: "NA" for i in column_names}

    writer = csv.DictWriter(sys.stdout, fieldnames=column_names, delimiter="\t")
    writer.writeheader()

    # perform analysis by group
    for group in group_to_truth_insertions:

        logger.info("-analyzing group: {}".format(group))

        truth_insertions = group_to_truth_insertions[group]
        pred_insertions = group_to_pred_insertions[group]

        analyze_insertions(group, truth_insertions, pred_insertions, writer, NA_row)

    sys.exit(0)


def extract_truth_breakpoint_info(row):

    # make chrA correspond to the ref genome
    if re.match("chr", row["chrA"]) and not re.match("chr", row["chrB"]):
        human_chr = row["chrA"]
        human_brkpt = (
            row["coordA_rend"] if row["orientA"] == "+" else row["coordA_lend"]
        )
        virus = row["chrB"]
        virus_brkpt = (
            row["coordB_lend"] if row["orientB"] == "+" else row["coordB_rend"]
        )

    elif re.match("chr", row["chrB"]) and not re.match("chr", row["chrA"]):
        human_chr = row["chrB"]
        human_brkpt = (
            row["coordB_lend"] if row["orientA"] == "+" else row["coordB_rend"]
        )
        virus = row["chrA"]
        virus_brkpt = (
            row["coordA_rend"] if row["orientB"] == "+" else row["coordA_lend"]
        )

    else:
        raise RuntimeError(
            "Error, cannot resolve ref genome chr from entry: {}".format(row)
        )

    return (human_chr, human_brkpt, virus, virus_brkpt)


class Insertion:
    def __init__(self, human_chr, human_brkpt, virus, virus_brkpt, row):

        self.name = "~".join([human_chr, str(human_brkpt), virus, str(virus_brkpt)])
        self.human_chr = str(human_chr)
        self.human_brkpt = int(human_brkpt)
        self.virus = str(virus)
        self.virus_brkpt = int(virus_brkpt)
        self.row = row


def analyze_insertions(group, truth_insertions, pred_insertions, writer, NA_row):

    truth_insertion_to_closest_preds = defaultdict(list)

    # map each pred insertion to the closest truth insertion and generate report
    for pred_insertion in pred_insertions:
        closest_dist = None
        closest_insertion = None

        for truth_insertion in truth_insertions:

            distance = None

            if (
                pred_insertion.virus == truth_insertion.virus
                and pred_insertion.human_chr == truth_insertion.human_chr
            ):

                distance = (
                    pred_insertion.virus_brkpt - truth_insertion.virus_brkpt
                ) ** 2 + abs(pred_insertion.human_brkpt - truth_insertion.human_brkpt)

                if closest_insertion is None or closest_dist > distance:
                    closest_insertion = truth_insertion
                    closest_dist = distance

        if closest_insertion is not None:
            truth_insertion_to_closest_preds[closest_insertion.name].append(
                pred_insertion
            )
            pred_insertion_dict = pred_insertion.row.to_dict()
            pred_insertion_dict["pred_insertion_name"] = pred_insertion.name
            pred_insertion_dict["truth_insertion_name"] = closest_insertion.name
            pred_insertion_dict["truth_human_chr"] = closest_insertion.human_chr
            pred_insertion_dict["truth_human_brkpt"] = closest_insertion.human_brkpt
            pred_insertion_dict["truth_virus"] = closest_insertion.virus
            pred_insertion_dict["truth_virus_brkpt"] = closest_insertion.virus_brkpt
            writer.writerow(pred_insertion_dict)

        else:
            pred_insertion_dict = pred_insertion.row.to_dict()
            pred_insertion_dict["pred_insertion_name"] = pred_insertion.name
            pred_insertion_dict["truth_insertion_name"] = "NA"
            pred_insertion_dict["truth_human_chr"] = "NA"
            pred_insertion_dict["truth_human_brkpt"] = "NA"
            pred_insertion_dict["truth_virus"] = "NA"
            pred_insertion_dict["truth_virus_brkpt"] = "NA"
            writer.writerow(pred_insertion_dict)

    # report remaining insertions that lack predictions assigned.
    for truth_insertion in truth_insertions:
        if truth_insertion.name not in truth_insertion_to_closest_preds:
            # FNs
            insertion_dict = NA_row
            insertion_dict["pred_insertion_name"] = "NA"
            insertion_dict["truth_insertion_name"] = truth_insertion.name
            insertion_dict["truth_human_chr"] = truth_insertion.row["human_chr"]
            insertion_dict["truth_human_brkpt"] = (truth_insertion.row["human_brkpt"],)
            insertion_dict["truth_virus"] = truth_insertion.row["virus"]
            insertion_dict["truth_virus_brkpt"] = truth_insertion.row["virus_brkpt"]
            writer.writerow(insertion_dict)

    return


if __name__ == "__main__":
    main()
