#!/usr/bin/env python3
# encoding: utf-8

import os, sys, re
import logging
import argparse
import pandas as pd
import pyranges as pr

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
        default=2,
        help="min sample count to define hotspot",
    )

    parser.add_argument(
        "--min_ev_reads",
        type=int,
        default=5,
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

    args = parser.parse_args()
    # args, unknown_args = parser.parse_known_args()

    if args.debug:
        logger.setLevel(logging.DEBUG)

    insertions_tsv_filename = args.insertions_tsv
    min_samples_per_hotspot = args.min_hotspot_samples
    window_size = args.window_size
    output_filename = args.output
    min_ev_reads = args.min_ev_reads

    data = pd.read_csv(insertions_tsv_filename, sep="\t")

    data = data[data["total"] >= min_ev_reads]

    columns = list(data.columns)
    if "#sample" in columns:
        data = data.rename(columns={"#sample": "sample"})
    elif "sample_name" in columns:
        data = data.rename(columns={"sample_name": "sample"})

    if "humanchr" in columns:
        data["human_chrom"] = data["humanchr"]
        data["virus_genome"] = data["virus"]
    else:
        data = data.apply(define_virus_and_genome_info, axis=1)

    hotspots = data.groupby("human_chrom").apply(find_hotspots, window_size=window_size)

    # join by proximity to defined hotspot

    data["Chromosome"] = data["human_chrom"]
    data["Start"] = data["human_coord"]
    data["End"] = data["human_coord"] + 1

    data = data.astype({"Start": "int", "End": "int"})

    data_pr = pr.PyRanges(data)
    hotspots_pr = pr.PyRanges(hotspots)

    data_nearest_pr = data_pr.nearest(hotspots_pr)

    data_w_hotspots = data_nearest_pr.df.copy()

    if min_samples_per_hotspot > 1:
        data_w_hotspots = data_w_hotspots[
            data_w_hotspots["hotspot_sample_counts"] >= min_samples_per_hotspot
        ]

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
            while (
                len(windowed_rows) > 0
                and human_coord - windowed_rows[0]["human_coord"] > window_size
            ):
                windowed_rows.pop(0)

        windowed_rows.append(row)
        uniq_samples = set()
        for x in windowed_rows:
            uniq_samples.add(x["sample"])

        sample_tallies.append([human_coord, len(uniq_samples)])

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


####################

if __name__ == "__main__":
    main()
