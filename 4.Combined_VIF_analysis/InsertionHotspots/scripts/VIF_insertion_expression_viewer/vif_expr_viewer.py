#!/usr/bin/env python3

import sys, os, re
import genomeview
from genomeview.intervaltrack import *
import argparse
import json
import bsddb3
import pysam
from collections import defaultdict
import seaborn as sns
import pyranges as pr
import csv
import subprocess

UTILDIR = os.path.join(os.path.dirname(__file__), "util")


def main():

    parser = argparse.ArgumentParser(
        description="generate expr view",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "--ref_annot_bed",
        type=str,
        required=True,
        help="reference gene structure annotation bed file",
    )
    parser.add_argument(
        "--ref_gene_spans_bed_tabix_gz",
        type=str,
        required=True,
        help="ref gene spans bed file bgzip'd and tabix indexed",
    )
    parser.add_argument(
        "--expr_matrix_bdbs",
        type=str,
        required=True,
        help="gene expr matrix stored in pickle format",
        nargs="+",
    )
    parser.add_argument(
        "--region",
        type=str,
        required=True,
        help="targeted region in igv/ucsc format:  'chr:lend-rend' ",
    )
    parser.add_argument(
        "--vif_insertions_tsv_tabix_gz",
        type=str,
        required=True,
        help="virus insertions tsv bgzip'd and tabix indexed",
    )
    parser.add_argument(
        "--output_filename",
        type=str,
        required=True,
        help="name of output file (can include pdf or svg extension",
    )
    parser.add_argument(
        "--gistic_matrix_bdbs",
        type=str,
        required=False,
        help="gistic gene-based CNV info in pickle form",
        nargs="+",
    )
    parser.add_argument(
        "--cnv_regions_tsv_tabix_gz",
        type=str,
        required=False,
        help="cnv regions in tsv bgzip'd and tabix indexed",
    )
    parser.add_argument(
        "--insertion_view_prefix",
        type=str,
        required=False,
        help="prefix for separate cnv insertion view pdf, included when --cnv_regions_tsv_tabix_gz is provided. Default: vif.{region}",
    )

    args = parser.parse_args()

    ref_annot_bed_filename = args.ref_annot_bed
    ref_gene_spans_bed = args.ref_gene_spans_bed_tabix_gz
    expr_matrix_bdbs = args.expr_matrix_bdbs
    region = args.region
    vif_insertions_tsv = args.vif_insertions_tsv_tabix_gz
    output_filename = args.output_filename
    # gistic_matrix_bdbs = args.gistic_matrix_bdbs
    cnv_regions_tsv = args.cnv_regions_tsv_tabix_gz
    insertion_view_prefix = args.insertion_view_prefix
    if insertion_view_prefix is None:
        insertion_view_prefix = f"vif.{region}"

    chrom, lend, rend = None, None, None

    m = re.search("^(\S+):(\d+)-(\d+)$", region)
    if m:
        chrom, lend, rend = m.groups()
        lend = int(lend)
        rend = int(rend)
    else:
        exit(f"cannot decipher {region} as chr:lend-rend formatting")

    viewport = {"chrom": chrom, "lend": lend, "rend": rend}

    genes_want = set()
    gene_spans_tabix = pysam.TabixFile(ref_gene_spans_bed)
    gene_spans = list()
    for locus in gene_spans_tabix.fetch(chrom, lend, rend):
        # chr1    491225  493241  RP4-669L17.8^ENSG00000250575.1  -
        print(locus)
        locus = locus.split()
        gene_sym = locus[-2].split("^")[0]
        locus[-2] = gene_sym
        gene_spans.append(locus)
        genes_want.add(gene_sym)

    insertions_tabix = pysam.TabixFile(vif_insertions_tsv)
    insertions = list()
    samples_want = set()
    for locus in insertions_tabix.fetch(chrom, lend, rend):
        locus = locus.split()
        insertions.append(locus)
        print(locus)
        # chr16   88696794        88696795        +       HPV16   881     -       37      LUSC-98-A538-TP chr16~88696794~+~HPV16~881~-
        samples_want.add(locus[-2])

    # retrieve expression info from the bdbs.

    gistic_info = None
    # if gistic_matrix_pickle_file:
    #    with open(gistic_matrix_pickle_file, "rb") as fh:
    #        gistic_info = pickle.load(fh)

    sample_to_cnv_regions = defaultdict(list)

    gene_expr_info = build_gene_expr_info(samples_want, genes_want, expr_matrix_bdbs)

    """
    print("GENE EXPR INFO:")
    for key in gene_expr_info.keys():
        print("key: " + str(key))
        print(str(gene_expr_info[key]))

    sys.exit(0)
    """

    if cnv_regions_tsv:

        # cnv info here: https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/CNV_Pipeline/
        # value: The GDC further transforms these copy number values into segment mean values, which are equal to log2(copy-number/ 2). Diploid regions will have a segment mean of zero, amplified regions will have positive values, and deletions will have negative values.

        cnv_tsv_outfile = f"{insertion_view_prefix}.cnv_info.tsv"
        cnv_tsv_ofh = open(cnv_tsv_outfile, "wt")
        tab_writer = csv.DictWriter(
            cnv_tsv_ofh,
            fieldnames=["chrom", "lend", "rend", "sample_name", "CN", "copy_number"],
            delimiter="\t",
        )
        tab_writer.writeheader()

        cnv_regions_tabix = pysam.TabixFile(cnv_regions_tsv)
        for cnv_region in cnv_regions_tabix.fetch(chrom, lend, rend):
            (
                cnv_region_chrom,
                cnv_region_lend,
                cnv_region_rend,
                sample_name,
                value,
            ) = cnv_region.split("\t")

            cnv_region_lend = int(cnv_region_lend)
            cnv_region_rend = int(cnv_region_rend)

            value = float(value)

            copy_number = 2 * 2 ** value

            cnv_dict = {
                "chrom": cnv_region_chrom,
                "lend": cnv_region_lend,
                "rend": cnv_region_rend,
                "copy_number": float(copy_number),
                "CN": value,
                "sample_name": sample_name,
            }

            if (
                chrom == viewport["chrom"]
                and cnv_region_lend < viewport["rend"]
                and cnv_region_rend > viewport["lend"]
            ):

                sample_to_cnv_regions[sample_name].append(cnv_dict)
                tab_writer.writerow(cnv_dict)

        cnv_tsv_ofh.close()

        # write the genespans and insertions to files too
        cnv_gene_spans_outfile = f"{insertion_view_prefix}.gene_spans.tsv"
        cnv_gene_expr_outfile = f"{insertion_view_prefix}.expr.tsv"
        cnv_gene_expr_ofh = open(cnv_gene_expr_outfile, "wt")
        with open(cnv_gene_spans_outfile, "wt") as ofh:
            tab_writer = csv.DictWriter(
                ofh,
                fieldnames=["chrom", "lend", "rend", "gene", "strand"],
                delimiter="\t",
            )
            tab_writer.writeheader()

            expr_tab_writer = csv.DictWriter(
                cnv_gene_expr_ofh,
                fieldnames=["gene", "sample_name", "expr", "rank"],
                delimiter="\t",
            )
            expr_tab_writer.writeheader()

            for gene_span in gene_spans:
                gene_chrom, gene_lend, gene_rend, gene_name, gene_strand = gene_span
                tab_writer.writerow(
                    {
                        "chrom": gene_chrom,
                        "lend": gene_lend,
                        "rend": gene_rend,
                        "gene": gene_name,
                        "strand": gene_strand,
                    }
                )

                for sample in gene_expr_info:
                    if gene_name in gene_expr_info[sample]:
                        expr_tab_writer.writerow(
                            {
                                "gene": gene_name,
                                "sample_name": sample,
                                "expr": gene_expr_info[sample][gene_name]["expr"],
                                "rank": gene_expr_info[sample][gene_name]["rank"],
                            }
                        )

        cnv_gene_expr_ofh.close()

        # write the insertions
        cnv_insertions_outfile = f"{insertion_view_prefix}.insertions.tsv"
        with open(cnv_insertions_outfile, "wt") as ofh:
            tab_writer = csv.DictWriter(
                ofh,
                fieldnames=[
                    "chrom",
                    "lend",
                    "rend",
                    "strand",
                    "virus",
                    "virus_brkpt",
                    "virus_strand",
                    "read_support",
                    "sample_name",
                    "insertion_name",
                ],
                delimiter="\t",
            )
            tab_writer.writeheader()

            for locus in insertions:
                (
                    hg_chrom,
                    hg_lend,
                    hg_rend,
                    hg_strand,
                    virus,
                    virus_brkpt,
                    virus_strand,
                    read_support,
                    sample_name,
                    insertion_name,
                ) = locus
                tab_writer.writerow(
                    {
                        "chrom": hg_chrom,
                        "lend": hg_lend,
                        "rend": hg_rend,
                        "strand": hg_strand,
                        "virus": virus,
                        "virus_brkpt": virus_brkpt,
                        "virus_strand": virus_strand,
                        "read_support": read_support,
                        "sample_name": sample_name,
                        "insertion_name": insertion_name,
                    }
                )

        # make cnv insertion view
        cnv_insertion_view_script = os.path.join(
            UTILDIR, "plot_CNV_insertions_for_cohort.Rscript"
        )
        cmd = " ".join(
            [
                cnv_insertion_view_script,
                f"--cnv_regions {cnv_tsv_outfile}",
                f"--gene_spans {cnv_gene_spans_outfile}",
                f"--viral_insertions {cnv_insertions_outfile}",
                f"--output_pdf {insertion_view_prefix}.cnv_insertion.plot.pdf",
            ]
        )

        # print(cmd, file=sys.stderr)
        # subprocess.check_call(cmd, shell=True)

        # make the cnv wilcoxon plot
        cnv_wilcoxon_plot_script = os.path.join(
            UTILDIR, "plot_CNV_for_cohort_hotspot.gene_Wilcoxons.Rscript"
        )
        cmd = " ".join(
            [
                cnv_wilcoxon_plot_script,
                f"--cnv_regions {cnv_tsv_outfile}",
                f"--gene_spans {cnv_gene_spans_outfile}",
                f"--viral_insertions {cnv_insertions_outfile}",
                f"--output_pdf {insertion_view_prefix}.cnv_insertion.wilcoxon_plot.pdf",
            ]
        )

        print(cmd, file=sys.stderr)
        subprocess.run(cmd, shell=True)

        # make cnv expr ranking heatmap insertion view:
        expr_ranking_insertion_view_script = os.path.join(
            UTILDIR, "plot_EXPR_outliers_for_cohort_hotspot.heatmap.Rscript"
        )
        cmd = " ".join(
            [
                expr_ranking_insertion_view_script,
                f"--expr_info {cnv_gene_expr_outfile}",
                f"--gene_spans {cnv_gene_spans_outfile}",
                f"--viral_insertions {cnv_insertions_outfile}",
                f"--output_pdf {insertion_view_prefix}.expr_ranking_insertion.plot.pdf",
            ]
        )

        # print(cmd, file=sys.stderr)
        # subprocess.check_call(cmd, shell=True)

        # make the expr wilcoxon plot:
        expr_wilcoxon_plot_script = os.path.join(
            UTILDIR, "plot_EXPR_outliers_for_cohort_hotspot.gene_Wilcoxons.Rscript"
        )
        cmd = " ".join(
            [
                expr_wilcoxon_plot_script,
                f"--expr_info {cnv_gene_expr_outfile}",
                f"--gene_spans {cnv_gene_spans_outfile}",
                f"--viral_insertions {cnv_insertions_outfile}",
                f"--output_pdf {insertion_view_prefix}.expr_ranking_insertion.wilcoxon_plot.pdf",
            ]
        )

        print(cmd, file=sys.stderr)
        subprocess.run(cmd, shell=True)

        # make the combined cnv, expression, insertion plot:
        expr_cnv_insertion_plotting_script = os.path.join(
            UTILDIR, "plot_CNV_and_EXPR_outliers_for_cohort_hotspot.heatmap.Rscript"
        )
        cmd = " ".join(
            [
                expr_cnv_insertion_plotting_script,
                f"--cnv_regions {cnv_tsv_outfile}",
                f"--expr_info {cnv_gene_expr_outfile}",
                f"--gene_spans {cnv_gene_spans_outfile}",
                f"--viral_insertions {cnv_insertions_outfile}",
                f"--output_pdf {insertion_view_prefix}.both_cnv_n_expr.plot.pdf",
            ]
        )

        print(cmd, file=sys.stderr)
        subprocess.check_call(cmd, shell=True)

    #####################
    # make genomeview img

    width = 1300
    doc = genomeview.Document(width)

    view = genomeview.GenomeView(chrom, lend, rend, "+")
    doc.add_view(view)

    axis_track = genomeview.Axis("axis")
    view.add_track(axis_track)

    bed_track = genomeview.BEDTrack(
        ref_annot_bed_filename, name=os.path.basename(ref_annot_bed_filename)
    )
    bed_track.draw_locus_labels = False
    view.add_track(bed_track)

    sample_to_insertions = defaultdict(list)

    for insertion in insertions:
        #  ['chr8', '128126453', '128126454', '+', 'HPV18', '2659', '-', '65', 'CESC-C5-A3HE-TP', 'chr8~128126453~+~HPV18~2659~-']
        (
            chrom,
            lend,
            rend,
            orient,
            virus,
            virus_coord,
            virus_orient,
            virus_read_count,
            sample_name,
            insertion_token,
        ) = insertion
        sample_to_insertions[sample_name].append(
            {
                "chrom": chrom,
                "lend": int(lend),
                "rend": int(rend),
                "virus": virus,
                "virus_coord": int(virus_coord),
                "virus_orient": virus_orient,
                "virus_read_count": int(virus_read_count),
            }
        )

    for sample_name, insertions_list in sample_to_insertions.items():
        # print("-PROCESSING: [" + sample_name + "]")
        # sample_name = "-".join(
        #    ["TCGA", *sample_name.split("-")[1:]]
        # )  # expr val samples all start with TCGA instead of the sample type
        # if re.search("-TP$", sample_name):
        #    sample_name = re.sub("-TP$", "-01A", sample_name)
        # elif re.search("-NT$", sample_name):
        #
        #    sample_name_test = re.sub("-NT$", "-01B", sample_name)
        #    if sample_name_test in gene_expr_info:
        #        sample_name = sample_name_test
        #    else:
        #        sample_name = re.sub("-NT$", "-11A", sample_name)

        print(f"******** sample_name: {sample_name}")

        expr_data = gene_expr_info[sample_name]

        print("]\texpr data for sample: {}".format(expr_data))

        gistic_vals = None
        if gistic_info:
            gistic_vals = gistic_info[sample_name]

        cnv_regions = None
        if sample_name in sample_to_cnv_regions:
            cnv_regions = sample_to_cnv_regions[sample_name]
            print("CNV regions for sample {}: {}".format(sample_name, cnv_regions))
        else:
            print("CNV - missing sample info for {}".format(sample_name))

        # print(str(expr_data))
        # exit("done")
        vif_sample_track = VIF_sample_track(
            sample_name,
            gene_spans,
            expr_data,
            insertions_list,
            gistic_vals,
            cnv_regions,
            viewport,
        )
        view.add_track(vif_sample_track)

    genomeview.save(doc, output_filename)

    sys.exit(0)


def build_gene_expr_info(samples_want, genes_want, expr_matrix_bdbs):

    gene_expr_info = defaultdict(dict)

    for expr_matrix_bdb in expr_matrix_bdbs:
        print("-exploring bdb for expr loading: {}".format(expr_matrix_bdb))
        dbase = bsddb3.hashopen(expr_matrix_bdb, "r")
        all_samples = set(json.loads(dbase[b"sample_list"]))
        all_genes = set(json.loads(dbase[b"gene_list"]))

        # print("all samples: {}".format(all_samples))

        # print("all genes: {}".format(all_genes))

        # get gene symbol mapping.
        symbol_to_gene_id = dict()
        for gene in all_genes:
            symbol = gene.split("^")[0]
            symbol_to_gene_id[symbol] = gene

        for sample in all_samples:
            for gene in genes_want:
                full_gene_id = symbol_to_gene_id[gene]
                # print("gene want: " + gene)
                if full_gene_id in all_genes:
                    sample_gene_key = sample + "::" + full_gene_id
                    expr_info = json.loads(dbase[sample_gene_key.encode()])
                    gene_expr_info[sample][gene] = expr_info

    return gene_expr_info


class VIF_sample_track(IntervalTrack):
    def __init__(
        self,
        sample_name,
        gene_intervals,
        expr_vals,
        insertion_positions,
        gistic_vals,
        cnv_regions,
        viewport,
    ):

        super().__init__([], name=sample_name)

        self.gene_intervals = gene_intervals
        self.expr_vals = expr_vals
        self.insertion_positions = sorted(
            insertion_positions, key=lambda x: x["virus_read_count"], reverse=True
        )
        self.max_virus_read_count = self.insertion_positions[0]["virus_read_count"]
        self.insertion_positions = sorted(insertion_positions, key=lambda x: x["lend"])

        self.gistic_vals = gistic_vals
        self.cnv_regions = cnv_regions

        self.viewport = viewport

        self.intervals = self

        self.draw_locus_labels = True
        self.include_locus_fn = None

        self.row_height = 12
        self.thick_width = self.row_height
        self.thin_width = 5
        self.thinnest_width = 1

        self.min_exon_width = 1

        self.num_colors = 10
        self.virus_colors = sns.color_palette("Reds", n_colors=self.num_colors).as_hex()

        self.max_gene_expr_rank = 50
        self.gene_expr_colors = sns.color_palette(
            "Blues", n_colors=self.num_colors
        ).as_hex()

    def __iter__(self):

        for gene_span in self.gene_intervals:
            # print("ladeda: " + str(gene_span))
            chrom, lend, rend, gene_name, orient = gene_span
            lend = int(lend)
            rend = int(rend)
            id = gene_name
            gene_symbol = gene_name.split("^")[0]
            interval = VIF_interval(
                id,
                chrom,
                lend,
                rend,
                orient,
                label=gene_symbol,
                type="gene",
                magnitude=None,
            )
            print("Interval: " + str(interval))
            yield interval

        seen = set()
        for virus_insertion in self.insertion_positions:
            print("virus: " + str(virus_insertion))
            (chrom, lend, rend, virus, virus_coord, virus_orient, virus_read_count) = (
                virus_insertion["chrom"],
                virus_insertion["lend"],
                virus_insertion["rend"],
                virus_insertion["virus"],
                virus_insertion["virus_coord"],
                virus_insertion["virus_orient"],
                virus_insertion["virus_read_count"],
            )

            id = ",".join([str(x) for x in [self.name, virus_insertion.values()]])
            if id in seen:
                continue
            interval = VIF_interval(
                id,
                chrom,
                int(lend),
                int(rend),
                virus_orient,
                label=virus + f"(C={virus_coord},T={virus_read_count})",
                type="virus",
                magnitude=int(virus_read_count),
            )
            yield interval
            seen.add(id)

        if self.cnv_regions is not None:
            for cnv_region in self.cnv_regions:
                print("cnv_region: " + str(cnv_region))
                (chrom, lend, rend, copy_number) = (
                    cnv_region["chrom"],
                    cnv_region["lend"],
                    cnv_region["rend"],
                    cnv_region["copy_number"],
                )
                id = f"{lend}-{rend}:{copy_number}"
                interval = VIF_interval(
                    id,
                    chrom,
                    int(lend),
                    int(rend),
                    None,
                    label=f"{copy_number:.2f}",
                    type="cnv",
                    magnitude=copy_number,
                )
                yield interval

        return

    def draw_gene(self, renderer, interval):
        # yield from super().draw_interval(renderer, interval)
        start = self.scale.topixels(interval.start)
        end = self.scale.topixels(interval.end)

        row = self.intervals_to_rows[interval.id]
        top = row * (self.row_height + self.margin_y)

        if interval.id in self.expr_vals:
            gene_expr_info = self.expr_vals[interval.id]
            rank_val = min(gene_expr_info["rank"], self.max_gene_expr_rank)
            color_idx = int(
                (self.max_gene_expr_rank - rank_val)
                / self.max_gene_expr_rank
                * (self.num_colors - 1)
            )
            print(f"color_idx: {color_idx}")
            color = self.gene_expr_colors[color_idx]
        else:
            print(f"-missing expr value for sample[{self.name}], gene[{interval.id}]")
            color = "#90EE90"  # light green for missing ones

        temp_label = interval.label
        if interval.label is None:
            temp_label = interval.id

        misc_opts = {"stroke": "none", "id": temp_label}

        if self.gistic_vals and interval.id in self.gistic_vals:
            print(f"got gistic for {interval.id} : {self.gistic_vals[interval.id]}")
            gene_gistic_info = self.gistic_vals[interval.id]
            gene_copy_number = int(gene_gistic_info["expr"])
            if gene_copy_number != 0:
                misc_opts["stroke"] = "red" if gene_copy_number > 0 else "blue"
                misc_opts["stroke-width"] = 2 * abs(gene_copy_number)

        if interval.strand is None:
            yield from renderer.rect(
                start, top, end - start, self.row_height, fill=color, **misc_opts
            )  # **{"stroke":"none", "id":temp_label})
        else:
            arrow_width = max(self.row_height, self.scale.relpixels(30))
            direction = "right" if interval.strand == "+" else "left"

            yield from renderer.block_arrow(
                start,
                top,
                end - start,
                self.row_height,
                arrow_width=arrow_width,
                direction=direction,
                fill=color,
                **misc_opts,
            )  # **{"stroke":"none", "id":temp_label})

        if interval.label is not None:
            yield from renderer.text(
                end + self.label_distance,
                top + self.row_height - 2,
                interval.label,
                anchor="start",
            )

    def draw_interval(self, renderer, interval):

        if interval.type == "gene":
            yield from self.draw_gene(renderer, interval)

        elif interval.type == "virus":
            yield from self.draw_virus(renderer, interval)

        elif interval.type == "cnv":
            yield from self.draw_cnv(renderer, interval)

        return

    def draw_cnv(self, renderer, interval):

        start = self.scale.topixels(interval.start)
        end = self.scale.topixels(interval.end) + 10

        magnitude = interval.magnitude

        print(f"pixels: {start}-{end}, magnitude {magnitude}, CNV region")

        row = self.intervals_to_rows[interval.id]
        top = row * (self.row_height + self.margin_y)

        # color = self.color_fn(interval)
        # color = self.virus_colors[ int ((magnitude / self.max_virus_read_count) * (self.num_colors - 1) ) ]
        color = "pink"

        temp_label = interval.label
        if interval.label is None:
            temp_label = interval.id

        yield from renderer.rect(
            start,
            top,
            end - start,
            self.row_height,
            fill=color,
            **{"stroke": "none", "id": temp_label},
        )

        if interval.label is not None:
            # yield from renderer.text(end+self.label_distance, top+self.row_height-2, interval.label, anchor="start")
            # place label based on viewport
            cnv_view_start = max(interval.start, self.viewport["lend"])
            cnv_view_end = min(interval.end, self.viewport["rend"])
            cnv_view_midpt = int((cnv_view_start + cnv_view_end) / 2)

            print(
                f"viewport [{self.viewport['lend']}-{self.viewport['rend']}] cnv_range[{interval.start}-{interval.end}] cnv midpt = {cnv_view_midpt}"
            )

            label_start_pixel_pos = self.scale.topixels(cnv_view_midpt)
            yield from renderer.text(
                label_start_pixel_pos,
                top + self.row_height - 2,
                interval.label,
                anchor="start",
            )

        return

    def draw_virus(self, renderer, interval):

        start = self.scale.topixels(interval.start)
        end = self.scale.topixels(interval.end) + 10

        magnitude = interval.magnitude

        print(
            f"pixels: {start}-{end}, magnitude {magnitude}, max read count {self.max_virus_read_count}"
        )

        row = self.intervals_to_rows[interval.id]
        top = row * (self.row_height + self.margin_y)

        # color = self.color_fn(interval)
        color = self.virus_colors[
            int((magnitude / self.max_virus_read_count) * (self.num_colors - 1))
        ]

        temp_label = interval.label
        if interval.label is None:
            temp_label = interval.id

        # yield from renderer.rect(start, top, end-start, self.row_height, fill=color,
        #     **{"stroke":"none", "id":temp_label})

        if interval.strand is None:
            # print("block")
            yield from renderer.rect(
                start,
                top,
                end - start,
                self.row_height,
                fill=color,
                **{"stroke": "none", "id": temp_label},
            )
        else:
            # print("arrow")
            arrow_width = max(
                self.row_height / 2, 0.1 * self.scale.relpixels(end - start)
            )
            print(f"arrow width: {arrow_width}")
            direction = "right" if interval.strand == "+" else "left"

            yield from renderer.block_arrow(
                start,
                top,
                end - start,
                self.row_height,
                arrow_width=arrow_width,
                direction=direction,
                fill=color,
                **{"stroke": "none", "id": temp_label},
            )

        if interval.label is not None:
            yield from renderer.text(
                end + self.label_distance,
                top + self.row_height - 2,
                interval.label,
                anchor="start",
            )


class VIF_interval(Interval):
    def __init__(
        self,
        id_,
        chrom,
        start,
        end,
        strand="+",
        label=None,
        type=["gene", "virus"][0],
        magnitude=None,
    ):

        super().__init__(id_, chrom, start, end, strand, label)

        self.type = type
        self.magnitude = magnitude


if __name__ == "__main__":
    main()
