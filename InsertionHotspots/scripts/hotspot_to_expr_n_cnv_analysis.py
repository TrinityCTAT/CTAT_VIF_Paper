#!/usr/bin/env python

import sys, os, re
import csv
from collections import defaultdict
import subprocess
import argparse
import bsddb3
import json
import pandas as pd

max_num_other_samples = 500

utildir = os.path.join(os.path.dirname(__file__), "util")

vif_expr_viewer_basedir = "~/GITHUB/CTAT_VIF/VIF_insertion_expression_viewer"


def main():

    parser = argparse.ArgumentParser(
        description="hotspot to expr and cnv analysis",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "--hotspots",
        type=str,
        required=True,
        help="hotspots file ie. hotspot_sample_data.threshold=5.for_expr.tsv",
    )
    parser.add_argument(
        "--ref_annot_bed",
        type=str,
        required=True,
        help="ref annotation gene structures in bed format",
    )
    parser.add_argument(
        "--ref_gene_spans", type=str, required=True, help="ref annot gene spans file"
    )
    parser.add_argument(
        "--expr_matrix",
        type=str,
        required=True,
        help="expression matrix(s) for cohort(s) data (eg. TCGA-CESC.htseq_fpkm-uq.genesym.tsv)",
        nargs="+",
    )
    parser.add_argument(
        "--min_hotspot_samples",
        type=int,
        default=5,
        help="min samples required per hotspot",
    )
    parser.add_argument(
        "--cnv_tsv", type=str, required=False, help="cnv tsv file(s)", nargs="*"
    )

    args = parser.parse_args()

    hotspots_filename = args.hotspots
    ref_annot_bed_filename = args.ref_annot_bed
    ref_gene_spans_filename = args.ref_gene_spans
    expr_matrix_filenames = args.expr_matrix
    min_hotspot_samples = args.min_hotspot_samples
    cnv_tsvs = args.cnv_tsv

    ## prep for hotspot genomeview
    hotspots_bed_gz = ensure_hotspots_bed_gz(hotspots_filename)
    ref_annot_bed_gz_filename = ensure_ref_annot_bed_gz_filename(ref_annot_bed_filename)
    ref_gene_spans_bed_gz = ensure_gene_spans_bed_gz(ref_gene_spans_filename)

    expr_matrix_bdbs = [ensure_expr_matrix_bdb(x) for x in expr_matrix_filenames]
    sample_core_to_samples_list, sample_to_bdb_dict = associate_bdb_with_samples(
        expr_matrix_bdbs
    )

    cnv_tsv_bed_gz = None
    if cnv_tsvs is not None:
        cnv_tsv_bed_gz = ensure_cnv_tabix(cnv_tsvs)

    print("-parsing gene spans", file=sys.stderr)
    gene_sym_to_genes = parse_gene_annots(ref_gene_spans_filename)

    all_hotspot_data = pd.read_csv(hotspots_filename, sep="\t")

    hotspot_to_samples = defaultdict(set)
    hotspot_to_genes = defaultdict(set)
    hotspot_to_gene_coords = defaultdict(list)

    with open(hotspots_filename) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:

            hotspot_tok = row["gene_regrouped_hotspot_name"]  # row["hotspot"]

            samplename = row["sample_id"]
            hotspot_to_samples[hotspot_tok].add(samplename)

            left_genes = row["left_genes"].split(";")
            right_genes = row["right_genes"].split(";")

            all_genes = left_genes + right_genes
            for gene in all_genes:
                if not re.search("RNA", gene):
                    hotspot_to_genes[hotspot_tok].add(gene)
                    gene_structs = gene_sym_to_genes[gene]
                    if len(gene_structs) == 1:
                        gene_struct = gene_structs[0]
                        if gene_struct["chrom"] == row["human_chrom"]:
                            hotspot_to_gene_coords[hotspot_tok].extend(
                                [gene_struct["lend"], gene_struct["rend"]]
                            )

    for hotspot in hotspot_to_samples:

        hotspot_samples = list(hotspot_to_samples[hotspot])

        if len(hotspot_samples) < min_hotspot_samples:
            continue

        hotspot_genes = list(hotspot_to_genes[hotspot])

        hotspot_chrom, hotspot_coord = hotspot.split("^")[0].split(":")

        if hotspot_chrom == "chrM":
            continue

        num_samples = len(hotspot_samples)
        num_genes = len(hotspot_genes)

        # skip those already processed for resume mode
        if os.path.exists(f"{hotspot}.{num_samples}s.expr_insertions_gview.pdf"):
            print(f"-skipping {hotspot}, already processed")
            continue

        # hotspot_data = all_hotspot_data[all_hotspot_data["hotspot"] == hotspot]
        hotspot_data = all_hotspot_data[
            all_hotspot_data["gene_regrouped_hotspot_name"] == hotspot
        ]
        hotspot_data.to_csv(
            f"{hotspot}.{num_samples}.hotspot_data.tsv", sep="\t", index=False
        )

        print(
            "\t".join(
                [
                    hotspot,
                    "sample_cnt: {}".format(num_samples),
                    "gene_cnt: {}".format(num_genes),
                    ",".join(hotspot_samples),
                    ",".join(hotspot_genes),
                ]
            )
        )

        hotspot_file_token = hotspot
        hotspot_file_token = hotspot_file_token.replace(";", "_")

        hotspot_expr_matrix_filename = (
            hotspot_file_token + f".{num_samples}s.expr.matrix"
        )
        hotspot_expr_matrix_samples_file = (
            hotspot_file_token + f".{num_samples}s.expr.samples.txt"
        )
        write_hotspot_expr_matrix(
            hotspot,
            hotspot_samples,
            hotspot_genes,
            gene_sym_to_genes,
            sample_core_to_samples_list,
            sample_to_bdb_dict,
            hotspot_expr_matrix_filename,
            hotspot_expr_matrix_samples_file,
        )

        actual_insertion_coords = hotspot_data["human_coord"].tolist()
        print("actual_insertion_coords: {}".format(actual_insertion_coords))
        hotspot_gene_coords = hotspot_to_gene_coords[hotspot]
        hotspot_gene_coords.extend(actual_insertion_coords)

        hotspot_gene_coords = sorted(hotspot_gene_coords)
        print("hotspot_gene_coords: {}".format(hotspot_gene_coords))

        flank = 1000
        hotspot_lend, hotspot_rend = (
            hotspot_gene_coords[0] - flank,
            hotspot_gene_coords[-1] + flank,
        )

        region_token = (
            #  f"{hotspot}.{num_samples}s.
            f"{hotspot_chrom}:{hotspot_lend}-{hotspot_rend}"
        )
        cmd = " ".join(
            [
                os.path.join(vif_expr_viewer_basedir, "vif_expr_viewer.py"),
                "--ref_annot_bed {}".format(ref_annot_bed_gz_filename),
                "--ref_gene_spans_bed {}".format(ref_gene_spans_bed_gz),
                "--expr_matrix_bdbs {}".format(" ".join(expr_matrix_bdbs)),
                "--region {}".format(region_token),
                "--vif_insertions_tsv_tabix_gz {}".format(hotspots_bed_gz),
                f"--insertion_view_prefix {hotspot_file_token}.{num_samples}s",
                f"--output_filename {hotspot_file_token}.{num_samples}s.expr_insertions_gview.pdf",
            ]
        )

        if cnv_tsv_bed_gz:
            cmd += " --cnv_regions_tsv_tabix_gz {}".format(cnv_tsv_bed_gz)

        # try:
        execute_cmd(cmd)

        # except:
        #    pass

        # continue

        # make heatmap
        # cmd = str("~/GITHUB/trinityrnaseq/Analysis/DifferentialExpression/PtR " +
        #          " -m {} -s {} ".format(expr_matrix_filename, expr_matrix_samples_file) +
        #          " --heatmap_scale_limits '-4,4' " +
        #          " --gene_clust none --heatmap --center_rows --order_columns_by_samples_file")

        # print("CMD: " + cmd)
        # subprocess.check_call(cmd, shell=True)

        # make expr ranking plot:
        # cmd = str(os.path.join(utildir,  "plot_expression_rankings.stacked_barplots.Rscript") +
        #          " --matrix {}".format(expr_matrix_filename) +
        #          " --samples {}".format(expr_matrix_samples_file) +
        #          " --outpng {}.expr_rank.stacked_barplots.png".format(hotspot) )

        # execute_cmd(cmd)

        # make expr ranking plot:
        cmd = str(
            os.path.join(utildir, "plot_expression_rankings.via_heatmap.Rscript")
            + " --matrix {}".format(hotspot_expr_matrix_filename)
            + " --samples {}".format(hotspot_expr_matrix_samples_file)
            + f" --outpng {hotspot_file_token}.{num_samples}s.expr_rank.heatmap.png"
        )

        execute_cmd(cmd)

    sys.exit(0)


def execute_cmd(cmd):
    print("CMD: " + cmd)
    subprocess.check_call(cmd, shell=True)


def parse_gene_annots(gene_spans_file):

    gene_sym_to_genes = defaultdict(list)

    with open(gene_spans_file) as fh:
        for line in fh:
            line = line.rstrip()
            ensg_id, chrom, lend, rend, orient, genesym, genetype = line.split("\t")
            lend = int(lend)
            rend = int(rend)
            midpt = int((lend + rend) / 2)

            gene_struct = {
                "genetok": genesym + "^" + ensg_id,
                "midpt": midpt,
                "chrom": chrom,
                "lend": lend,
                "rend": rend,
            }

            gene_sym_to_genes[genesym].append(gene_struct)

    return gene_sym_to_genes


def write_hotspot_expr_matrix(
    hotspot,
    hotspot_samples,
    hotspot_genes,
    gene_sym_to_genes,
    sample_core_to_samples_list,
    sample_to_bdb_dict,
    expr_matrix_filename,
    expr_matrix_samples_file,
):

    print(
        "-writing hotspot expr matrix: {}".format(expr_matrix_filename), file=sys.stderr
    )

    chrom, hotspot_coord = hotspot.split("^")[0].split(":")
    hotspot_coord = int(hotspot_coord)

    samples_want = set()
    gene_toks_want = set()

    samples_ofh = open(expr_matrix_samples_file, "wt")

    hotspot_sample_cores = set()
    for hotspot_sample in hotspot_samples:
        hotspot_sample_core = "-".join(hotspot_sample.split("-")[1:3])
        hotspot_sample_cores.add(hotspot_sample_core)
        for samplename in sample_core_to_samples_list[hotspot_sample_core]:
            samples_want.add(samplename)
            print("\t".join(["HPVins", samplename]), file=samples_ofh)

    other_sample_cores = set(sample_core_to_samples_list.keys())
    other_sample_cores = list(other_sample_cores - set(hotspot_sample_cores))
    other_sample_cores = other_sample_cores[0:max_num_other_samples]
    for other_sample_core in other_sample_cores:
        for samplename in sample_core_to_samples_list[other_sample_core]:
            samples_want.add(samplename)
            print("\t".join(["other", samplename]), file=samples_ofh)

    samples_ofh.close()

    gene_tok_to_midpt = dict()

    for hotspot_gene in hotspot_genes:
        genes_list = gene_sym_to_genes[hotspot_gene]
        for gene in genes_list:
            if (
                gene["chrom"] == chrom
                and abs(gene["midpt"] - hotspot_coord) <= 10e6
                and gene["genetok"]
                not in [
                    "IGH@-ext^IGH.g@-ext",
                    "IGH-@-ext^IGH-.g@-ext",
                    "IGL@-ext_IGL-@-ext",
                    "IGL-@-ext^IGL-.g@-ext",
                    "IGL@-ext^IGL.g@-ext",
                ]
            ):

                gene_tok = gene["genetok"]
                gene_toks_want.add(gene_tok)
                gene_tok_to_midpt[gene_tok] = gene["midpt"]

    gene_toks_want = list(gene_toks_want)
    gene_toks_want = sorted(gene_toks_want, key=lambda x: gene_tok_to_midpt[x])

    samples_want = list(samples_want)
    # build matrix.
    with open(expr_matrix_filename, "wt") as ofh:
        print("\t" + "\t".join(samples_want), file=ofh)
        for gene_tok in gene_toks_want:
            vals = [gene_tok]
            found_all_expr_vals = True
            for sample_name in samples_want:
                expr_val = get_gene_expr_val(
                    gene_tok, sample_name, sample_to_bdb_dict
                )  # gene_tok_to_sample_expr[gene_tok][sample_name]
                if expr_val is not None:
                    vals.append(expr_val)
                else:
                    found_all_expr_vals = False
                    break

            if not found_all_expr_vals:
                print(
                    "-warning, gene_tok: {} not found in expr matrix.".format(gene_tok),
                    file=sys.stderr,
                )
                continue
            else:
                print("\t".join([str(x) for x in vals]), file=ofh)

    return


def ensure_cnv_tabix(cnv_tsv_filenames):

    assert type(cnv_tsv_filenames) == list, "Error, cnv_tsv_filenames must be a list"

    cnv_bed = "cnv_info.bedlike.tsv"
    cnv_bed_gz = cnv_bed + ".gz"
    cnv_bed_gz_tabix = cnv_bed_gz + ".tbi"

    if not os.path.exists(cnv_bed_gz) and not os.path.exists(cnv_bed_gz_tabix):

        cmd = " ".join(
            [
                os.path.join(vif_expr_viewer_basedir, "util/cnv_to_bedlike_tsv.py"),
                "--cnv_tsv",
                " ".join(cnv_tsv_filenames),
                "--output {}".format(cnv_bed),
            ]
        )

        execute_cmd(cmd)

        execute_cmd(f"bgzip {cnv_bed}")

        execute_cmd(f"tabix -p bed {cnv_bed_gz}")

    return cnv_bed_gz


####
def ensure_hotspots_bed_gz(hotspots_filename):

    hotspots_bed = hotspots_filename + ".bedlike.tsv"
    hotspots_bed_gz = hotspots_bed + ".gz"
    hotspots_bed_gz_tabix = hotspots_bed_gz + ".tbi"

    if not os.path.exists(hotspots_bed_gz) and not os.path.exists(
        hotspots_bed_gz_tabix
    ):

        cmd = " ".join(
            [
                os.path.join(
                    vif_expr_viewer_basedir, "util/hotspots_to_bedlike_tsv.py"
                ),
                hotspots_filename,
            ]
        )

        execute_cmd(cmd)

        execute_cmd(f"bgzip {hotspots_bed}")

        execute_cmd(f"tabix -p bed {hotspots_bed_gz}")

    return hotspots_bed_gz


####
def ensure_gene_spans_bed_gz(ref_gene_spans_filename):

    ref_gene_spans_bed_gz = ref_gene_spans_filename + ".bed.gz"

    if not os.path.exists(ref_gene_spans_bed_gz):
        cmd = " ".join(
            [
                os.path.join(vif_expr_viewer_basedir, "util/gtf_gene_spans_to_bed.py"),
                ref_gene_spans_filename,
                f" | sort -k 1,1 -k2,2g -k3,3g | bgzip -c > {ref_gene_spans_bed_gz}",
            ]
        )
        execute_cmd(cmd)

        execute_cmd(f"tabix {ref_gene_spans_bed_gz}")

    return ref_gene_spans_bed_gz


####
def ensure_ref_annot_bed_gz_filename(ref_annot_bed_filename):

    ref_annot_bed_gz_filename = ref_annot_bed_filename + ".sorted.bed.gz"
    ref_annot_bed_gz_tabix_filename = ref_annot_bed_gz_filename + ".tbi"

    if not os.path.exists(ref_annot_bed_gz_filename) and not os.path.exists(
        ref_annot_bed_gz_tabix_filename
    ):
        cmd = f"sort -k1,1 -k2,2g -k3,3g {ref_annot_bed_filename} | bgzip -c > {ref_annot_bed_gz_filename}"
        execute_cmd(cmd)

        execute_cmd(f"tabix {ref_annot_bed_gz_filename}")

    return ref_annot_bed_gz_filename


####
def ensure_expr_matrix_bdb(expr_matrix_filename):

    expr_matrix_bdb = expr_matrix_filename + ".bdb"

    if not os.path.exists(expr_matrix_bdb):
        cmd = " ".join(
            [
                os.path.join(
                    vif_expr_viewer_basedir, "util/index_expression_matrix.bsddb3.py"
                ),
                expr_matrix_filename,
            ]
        )
        execute_cmd(cmd)

    return expr_matrix_bdb


####
def associate_bdb_with_samples(expr_matrix_bdbs_list):

    sample_to_bdb = dict()
    sample_core_to_samples_list = defaultdict(list)

    for bdb in expr_matrix_bdbs_list:
        dbase = bsddb3.hashopen(bdb, "r")
        sample_list = json.loads(dbase[b"sample_list"])
        for sample in sample_list:
            sample_to_bdb[sample] = bdb
            sample_core = "-".join(sample.split("-")[1:3])
            sample_core_to_samples_list[sample_core].append(sample)

    return sample_core_to_samples_list, sample_to_bdb


####
def get_gene_expr_val(gene_tok, sample_name, sample_to_bdb_dict):

    bdb = sample_to_bdb_dict[sample_name]
    dbase = bsddb3.hashopen(bdb, "r")
    sample_gene_token = sample_name + "::" + gene_tok
    expr_info = json.loads(dbase[sample_gene_token.encode()])
    expr_val = expr_info["expr"]
    return expr_val


if __name__ == "__main__":
    main()
