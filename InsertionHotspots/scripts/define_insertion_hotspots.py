#!/usr/bin/env python3
# encoding: utf-8

import os, sys, re
import logging
import argparse
import pandas as pd
import pyranges as pr
import math
import statistics
import random

if sys.version_info[0] != 3:
    print("This script requires Python 3")
    exit(1)


logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s : %(levelname)s : %(message)s",
    datefmt="%H:%M:%S",
)
logger = logging.getLogger(__name__)


def main():

    parser = argparse.ArgumentParser(
        description="define hotspot insertions",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "--insertions_tsv", type=str, required=True, help="define insertions hotspots"
    )

    parser.add_argument(
        "--min_hotspot_samples",
        type=int,
        default=1,
        help="min sample count to define hotspot",
    )

    parser.add_argument(
        "--min_ev_reads",
        type=int,
        default=2,
        help="min number of evidence reads for each insertion candidate",
    )

    parser.add_argument(
        "--window_size",
        type=int,
        required=False,
        default=1000000,
        help="window in which to define hotspots",
    )
    parser.add_argument(
        "--debug", action="store_true", default=False, help="debug mode"
    )

    parser.add_argument(
        "--output",
        type=str,
        required=True,
        help="output with insertions assigned to hotspots",
    )

    parser.add_argument(
        "--randomize",
        action="store_true",
        default=False,
        help="randomize insertion positions",
    )

    args = parser.parse_args()
    # args, unknown_args = parser.parse_known_args()

    if args.debug:
        logger.setLevel(logging.DEBUG)

    insertions_tsv_filename = args.insertions_tsv
    min_samples_per_hotspot = args.min_hotspot_samples
    window_size = args.window_size
    output_filename = args.output
    min_ev_reads = args.min_ev_reads

    data = pd.read_csv(insertions_tsv_filename, sep="\t", low_memory=False)
    logger.info("Num insertions input: {}".format(data.shape[0]))

    data = data[data["total"] >= min_ev_reads]
    logger.info(
        "After filtering for min {} evidence read pairs have {} insertions".format(
            min_ev_reads, data.shape[0]
        )
    )

    columns = list(data.columns)
    if "#sample" in columns:
        data = data.rename(columns={"#sample": "sample"})
    elif "sample_name" in columns:
        data = data.rename(columns={"sample_name": "sample"})

    if "humanchr" in columns:
        data = data.rename(columns={"humanchr": "human_chrom", "virus": "virus_genome"})
    else:
        data = data.apply(define_virus_and_genome_info, axis=1)

    # distill initially to top insertion within window_size per sample.
    hotspot_target_data = data.copy()
    hotspot_target_data["insertion_locus_coord"] = round(
        hotspot_target_data["human_coord"] / window_size
    )

    hotspot_target_data = (
        hotspot_target_data.groupby(
            ["participant", "human_chrom", "insertion_locus_coord"]
        )
        .apply(lambda x: x.sort_values("total_rpm", ascending=False).head(1))
        .reset_index(drop=True)
    )

    logger.info(
        "After distilling to top insertion per participant in each window, have {} insertions to explore for hotspots".format(
            hotspot_target_data.shape[0]
        )
    )

    if args.randomize:
        logger.info("-Randomizing insertion positions")
        randomizer = Random_genome_position()
        hotspot_target_data = randomizer.randomize_insertion_positions(
            hotspot_target_data
        )
        data = hotspot_target_data.copy()

    logger.info("Finding hotspots")
    hotspots = hotspot_target_data.groupby("human_chrom").apply(
        find_hotspots, window_size=window_size
    )

    # join by proximity to defined hotspot

    data["Chromosome"] = data["human_chrom"]
    data["Start"] = data["human_coord"]
    data["End"] = data["human_coord"] + 1

    data = data.astype({"Start": "int", "End": "int"})

    data_pr = pr.PyRanges(data)
    hotspots_pr = pr.PyRanges(hotspots)

    data_nearest_pr = data_pr.nearest(hotspots_pr)

    data_w_hotspots = data_nearest_pr.df.copy()
    logger.info(
        "Number of insertions after assigning hotspots: {}".format(
            data_w_hotspots.shape[0]
        )
    )

    # reassign sample counts

    def redefine_hotspot_sample_counts(hotspot_df):
        hotspot_samples = sorted(list(hotspot_df["sample"].unique()))
        hotspot_df["hotspot_sample_counts_redef"] = len(hotspot_samples)
        hotspot_df["hotspot_sample_counts_redef_list"] = ",".join(hotspot_samples)
        return hotspot_df

    data_w_hotspots = data_w_hotspots.groupby("hotspot").apply(
        redefine_hotspot_sample_counts
    )

    if min_samples_per_hotspot > 1:
        logger.info(
            "Restricting to those with at least {} samples per hotspot".format(
                min_samples_per_hotspot
            )
        )
        data_w_hotspots = data_w_hotspots[
            data_w_hotspots["hotspot_sample_counts_redef"] >= min_samples_per_hotspot
        ]
        logger.info(
            "Number of insertions after filtering: {}".format(data_w_hotspots.shape[0])
        )

    data_w_hotspots.to_csv(output_filename, sep="\t", index=False)

    sys.exit(0)


def define_virus_and_genome_info(row):

    chrA, coordA = row["chrA"], row["coordA"]
    chrB, coordB = row["chrB"], row["coordB"]

    assert (re.match("^chr[\\dMXY]+", chrA) is not None) ^ (
        re.match("^chr[\\dMXY]+", chrB) is not None
    ), "Error, cannot figure out which of pair is human chromosome"

    if re.match("^chr[\\dMXY]+", chrA):
        human_chrom, human_coord = chrA, coordA
        virus_genome, virus_coord = chrB, coordB

    else:
        human_chrom, human_coord = chrB, coordB
        virus_genome, virus_coord = chrA, coordA

    row["human_chrom"] = human_chrom
    row["human_coord"] = human_coord
    row["virus_genome"] = virus_genome
    row["virus_coord"] = virus_coord

    return row


def find_hotspots(chrom_df, window_size):

    chrom_df = chrom_df.sort_values(["human_coord"])

    sample_tallies = []

    windowed_rows = []

    for i, row in chrom_df.iterrows():

        human_coord = row["human_coord"]

        if windowed_rows:
            # remove entries outside of window.
            while (
                len(windowed_rows) > 0
                and human_coord - windowed_rows[0]["human_coord"] > window_size
            ):
                windowed_rows.pop(0)

        windowed_rows.append(row)
        uniq_samples = set()
        for x in windowed_rows:
            uniq_samples.add(x["sample"])

        window_gene_midpt = round(
            statistics.mean([row["human_coord"] for row in windowed_rows])
        )
        sample_tallies.append([window_gene_midpt, len(uniq_samples)])

    # define the hotspots
    selected_hotspots = []

    sample_tallies = sorted(sample_tallies, key=lambda x: x[1], reverse=True)

    for hotspot_candidate in sample_tallies:
        coordinate, sample_count = hotspot_candidate
        if not within_range_of_hotspot(
            hotspot_candidate, selected_hotspots, window_size
        ):
            selected_hotspots.append(hotspot_candidate)

    hotspot_coordinates = [x[0] for x in selected_hotspots]
    hotspot_sample_counts = [x[1] for x in selected_hotspots]

    human_chrom = chrom_df["human_chrom"].unique()[0]

    hotspot_df = pd.DataFrame(
        {
            "Chromosome": [human_chrom] * len(hotspot_coordinates),
            "Start": hotspot_coordinates,
            "End": [x + 1 for x in hotspot_coordinates],
            "hotspot": [f"{human_chrom}:{x}" for x in hotspot_coordinates],
            "hotspot_sample_counts": hotspot_sample_counts,
        }
    )

    hotspot_df = hotspot_df.astype({"Start": "int", "End": "int"})

    return hotspot_df


def within_range_of_hotspot(hotspot_candidate, selected_hotspots, window_size):

    for selected_hotspot in selected_hotspots:
        if abs(hotspot_candidate[0] - selected_hotspot[0]) < window_size:
            return True

    return False


class Random_genome_position:
    def __init__(self):

        chr_lengths = {
            "chr1": 248956422,
            "chr2": 242193529,
            "chr3": 198295559,
            "chr4": 190214555,
            "chr5": 181538259,
            "chr6": 170805979,
            "chr7": 159345973,
            "chr8": 145138636,
            "chr9": 138394717,
            "chr10": 133797422,
            "chr11": 135086622,
            "chr12": 133275309,
            "chr13": 114364328,
            "chr14": 107043718,
            "chr15": 101991189,
            "chr16": 90338345,
            "chr17": 83257441,
            "chr18": 80373285,
            "chr19": 58617616,
            "chr20": 64444167,
            "chr21": 46709983,
            "chr22": 50818468,
            "chrX": 156040895,
            "chrY": 57227415,
        }

        chromosomes = list(chr_lengths.keys())

        self.sum_genome_length = sum(chr_lengths.values())

        self.ordered_chromosomes = list()
        cumsum = 0
        for chromosome in chromosomes:
            chr_len = chr_lengths[chromosome]
            self.ordered_chromosomes.append(
                {"chrom": chromosome, "lend": cumsum + 1, "rend": cumsum + chr_len}
            )
            cumsum += chr_len

    def randomize_insertion_positions(self, data):
        def assign_random_chr_position(row):

            rand_chrom_pos = random.randint(1, self.sum_genome_length)

            for chrom_struct in self.ordered_chromosomes:
                if (
                    rand_chrom_pos >= chrom_struct["lend"]
                    and rand_chrom_pos <= chrom_struct["rend"]
                ):
                    chrom = chrom_struct["chrom"]
                    chrom_pos = rand_chrom_pos - chrom_struct["lend"] + 1

                    row["human_chrom"] = chrom
                    row["human_coord"] = chrom_pos

                    return row

        data = data.apply(assign_random_chr_position, axis=1)

        return data


####################

if __name__ == "__main__":
    main()
