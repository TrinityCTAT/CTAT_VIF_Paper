#!/usr/bin/env python

import sys, os, re
import csv
from collections import defaultdict
import subprocess
import argparse

max_num_other_samples = 500

utildir = os.path.join(os.path.dirname(__file__), "util")

vif_expr_viewer_basedir = "~/GITHUB/CTAT_VIF/VIF_insertion_expression_viewer"

def main():


    parser = argparse.ArgumentParser(description="assign neighboring genes to hotspots", formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("--hotspots", type=str, required=True, help="hotspots file ie. hotspot_sample_data.threshold=5.for_expr.tsv")
    parser.add_argument("--ref_annot_bed", type=str, required=True, help="ref annotation gene structures in bed format")
    parser.add_argument("--ref_gene_spans", type=str, required=True, help='ref annot gene spans file')
    parser.add_argument("--expr_matrix", type=str, required=True, help='expression matrix for cohort data (eg. TCGA-CESC.htseq_fpkm-uq.genesym.tsv)')
    parser.add_argument("--min_hotspot_samples", type=int, default=5, help='min samples required per hotspot')
    parser.add_argument("--cnv_tsv", type=str, required=False, help="cnv tsv file")
    
    args = parser.parse_args()

    hotspots_filename = args.hotspots
    ref_annot_bed_filename = args.ref_annot_bed
    ref_gene_spans_filename = args.ref_gene_spans
    expr_matrix_filename= args.expr_matrix
    min_hotspot_samples = args.min_hotspot_samples
    cnv_tsv = args.cnv_tsv

    ## prep for hotspot genomeview
    hotspots_bed_gz = ensure_hotspots_bed_gz(hotspots_filename)
    ref_annot_bed_gz_filename = ensure_ref_annot_bed_gz_filename(ref_annot_bed_filename)
    ref_gene_spans_bed_gz = ensure_gene_spans_bed_gz(ref_gene_spans_filename)
    expr_matrix_pickle = ensure_expr_matrix_pickle(expr_matrix_filename)
    cnv_tsv_bed_gz = ensure_cnv_bed_gz(cnv_tsv) if cnv_tsv else None

    
    print("-parsing gene spans", file=sys.stderr)
    gene_sym_to_genes = parse_gene_annots(ref_gene_spans_filename)

    print("-parsing expr matrix", file=sys.stderr)
    sample_core_to_samples_list, gene_tok_to_sample_expr = parse_gene_expr_matrix(expr_matrix_filename)
    
    
    hotspot_to_samples = defaultdict(set)
    hotspot_to_genes = defaultdict(set)
    hotspot_to_gene_coords = defaultdict(list)
    
    with open(hotspots_filename) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            
            hotspot_tok = row['hotspot']
            
            samplename = row['sample']
            hotspot_to_samples[hotspot_tok].add(samplename)

            left_genes = row['left_genes'].split(";")
            right_genes = row['right_genes'].split(";")

            all_genes = left_genes + right_genes
            for gene in all_genes:
                if not re.search("RNA", gene):
                    hotspot_to_genes[hotspot_tok].add(gene)
                    gene_structs = gene_sym_to_genes[gene]
                    if len(gene_structs) == 1:
                        gene_struct = gene_structs[0]
                        if gene_struct['chrom'] == row['human_chrom']:
                            hotspot_to_gene_coords[hotspot_tok].extend([gene_struct['lend'], gene_struct['rend'] ])

    
    for hotspot in hotspot_to_samples:
        hotspot_samples = list(hotspot_to_samples[hotspot])

        if len(hotspot_samples) < min_hotspot_samples:
            continue
        
        hotspot_genes = list(hotspot_to_genes[hotspot])

        hotspot_chrom, hotspot_coord = hotspot.split(':')

        if hotspot_chrom == "chrM":
            continue
        
        hotspot_gene_coords = hotspot_to_gene_coords[hotspot]
        hotspot_gene_coords = sorted(hotspot_gene_coords)
        hotspot_lend, hotspot_rend = hotspot_gene_coords[0], hotspot_gene_coords[-1]
                
        num_samples = len(hotspot_samples)
        num_genes = len(hotspot_genes)

        print("\t".join([hotspot,
                         "sample_cnt: {}".format(num_samples),
                         "gene_cnt: {}".format(num_genes),
                         ",".join(hotspot_samples),
                         ",".join(hotspot_genes)]))
            
        expr_matrix_filename = hotspot + ".expr.matrix"
        expr_matrix_samples_file = hotspot + ".expr.samples.txt"
        write_hotspot_expr_matrix(hotspot, hotspot_samples, hotspot_genes,
                                  gene_sym_to_genes,
                                  sample_core_to_samples_list, gene_tok_to_sample_expr,
                                  expr_matrix_filename, expr_matrix_samples_file)

        
        region_token = f"{hotspot_chrom}:{hotspot_lend}-{hotspot_rend}"
        cmd = " ".join([os.path.join(vif_expr_viewer_basedir, "vif_expr_viewer.py"),
                        "--ref_annot_bed {}".format(ref_annot_bed_gz_filename),
                        "--ref_gene_spans_bed {}".format(ref_gene_spans_bed_gz),
                        "--expr_matrix_pickle {}".format(expr_matrix_pickle),
                        "--region {}".format(region_token),
                        "--vif_insertions_tsv_tabix_gz {}".format(hotspots_bed_gz),
                        "--output_filename {}.expr_insertions_gview.pdf".format(hotspot) ] )

        if cnv_tsv_bed_gz:
            cmd += f" --cnv_regions_tsv_tabix_gz {cnv_tsv_bed_gz}"

        #try:
        execute_cmd(cmd)

        #except:
        #    pass
        
        #continue
        
        # make heatmap
        #cmd = str("~/GITHUB/trinityrnaseq/Analysis/DifferentialExpression/PtR " +
        #          " -m {} -s {} ".format(expr_matrix_filename, expr_matrix_samples_file) +
        #          " --heatmap_scale_limits '-4,4' " +
        #          " --gene_clust none --heatmap --center_rows --order_columns_by_samples_file")
        
        #print("CMD: " + cmd)
        #subprocess.check_call(cmd, shell=True)

       
        # make expr ranking plot:
        #cmd = str(os.path.join(utildir,  "plot_expression_rankings.stacked_barplots.Rscript") +
        #          " --matrix {}".format(expr_matrix_filename) +
        #          " --samples {}".format(expr_matrix_samples_file) +
        #          " --outpng {}.expr_rank.stacked_barplots.png".format(hotspot) )

        #execute_cmd(cmd)
       
        
        # make expr ranking plot:
        cmd = str(os.path.join(utildir, "plot_expression_rankings.via_heatmap.Rscript") +
                  " --matrix {}".format(expr_matrix_filename) +
                  " --samples {}".format(expr_matrix_samples_file) +
                  " --outpng {}.expr_rank.heatmap.png".format(hotspot) )

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
            midpt = int( (lend+rend)/2)
            
            gene_struct = { 'genetok' : genesym + "^" + ensg_id,
                            'midpt' : midpt,
                            'chrom' : chrom,
                            'lend' : lend,
                            'rend' : rend }

            gene_sym_to_genes[genesym].append(gene_struct)

    return gene_sym_to_genes


def parse_gene_expr_matrix(expr_matrix):

    sample_core_to_samples_list = defaultdict(list)

    gene_tok_to_sample_expr = defaultdict(dict)
    
    
    with open(expr_matrix) as fh:
        sample_headers = next(fh)
        sample_headers = sample_headers.rstrip()
        samples = sample_headers.split("\t")
        samples.pop(0)
        
        for samplename in samples:
            sample_core = "-".join(samplename.split("-")[1:3])
            sample_core_to_samples_list[sample_core].append(samplename)

        
        for expr_row in fh:
            expr_row = expr_row.rstrip()
            expr_vals = expr_row.split("\t")
            gene_tok = expr_vals.pop(0)
            for i, expr_val in enumerate(expr_vals):
                samplename = samples[i]
                gene_tok_to_sample_expr[gene_tok][samplename] = expr_val

    return sample_core_to_samples_list, gene_tok_to_sample_expr



def write_hotspot_expr_matrix(hotspot, hotspot_samples, hotspot_genes,
                                  gene_sym_to_genes,
                                  sample_core_to_samples_list, gene_tok_to_sample_expr,
                                  expr_matrix_filename, expr_matrix_samples_file):

    print("-writing expr matrix: {}".format(expr_matrix_filename), file=sys.stderr)
    
    chrom, hotspot_coord = hotspot.split(":")
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
            if (gene['chrom'] == chrom and
                abs(gene['midpt'] - hotspot_coord) <=  10e6):
                gene_tok = gene['genetok']
                gene_toks_want.add(gene_tok)
                gene_tok_to_midpt[gene_tok] = gene['midpt']
                
    gene_toks_want = list(gene_toks_want)
    gene_toks_want = sorted(gene_toks_want, key=lambda x: gene_tok_to_midpt[x])

    samples_want = list(samples_want)
    # build matrix.
    with open(expr_matrix_filename, "wt") as ofh:
        print("\t" + "\t".join(samples_want), file=ofh)
        for gene_tok in gene_toks_want:
            if gene_tok not in gene_tok_to_sample_expr:
                print("-warning, gene_tok: {} not found in expr matrix.".format(gene_tok), file=sys.stderr)
                continue
            
            vals = [gene_tok]
            for sample_name in samples_want:
                expr_val = gene_tok_to_sample_expr[gene_tok][sample_name]
                vals.append(expr_val)

            print("\t".join(vals), file=ofh)

    return




def ensure_cnv_bed_gz(cnv_tsv_filename):

    cnv_bed = cnv_tsv_filename + ".bedlike.tsv"
    cnv_bed_gz = cnv_bed + ".gz"
    cnv_bed_gz_tabix = cnv_bed_gz + ".tbi"

    if not os.path.exists(cnv_bed_gz) and not os.path.exists(cnv_bed_gz_tabix):

        cmd = " ".join([os.path.join(vif_expr_viewer_basedir, "util/cnv_to_bedlike_tsv.py"),
                        cnv_tsv_filename])

        execute_cmd(cmd)

        execute_cmd(f"bgzip {cnv_bed}")

        execute_cmd(f"tabix -p bed {cnv_bed_gz}")

    
    return cnv_bed_gz    


####
def ensure_hotspots_bed_gz(hotspots_filename):

    hotspots_bed = hotspots_filename + ".bedlike.tsv"
    hotspots_bed_gz = hotspots_bed + ".gz"
    hotspots_bed_gz_tabix = hotspots_bed_gz + ".tbi"


    if not os.path.exists(hotspots_bed_gz) and not os.path.exists(hotspots_bed_gz_tabix):
    
        cmd = " ".join([os.path.join(vif_expr_viewer_basedir, "util/hotspots_to_bedlike_tsv.py"),
                    hotspots_filename])
        
        execute_cmd(cmd)

        execute_cmd(f"bgzip {hotspots_bed}")

        execute_cmd(f"tabix -p bed {hotspots_bed_gz}")


    return hotspots_bed_gz


####
def ensure_gene_spans_bed_gz(ref_gene_spans_filename):

    ref_gene_spans_bed_gz = ref_gene_spans_filename + ".bed.gz"

    if not os.path.exists(ref_gene_spans_bed_gz):
        cmd = " ".join([os.path.join(vif_expr_viewer_basedir, "util/gtf_gene_spans_to_bed.py"),
                        ref_gene_spans_filename,
                        f" | sort -k 1,1 -k2,2g -k3,3g | bgzip -c > {ref_gene_spans_bed_gz}"])
        execute_cmd(cmd)

        execute_cmd(f"tabix {ref_gene_spans_bed_gz}")


    return ref_gene_spans_bed_gz
        
####
def ensure_ref_annot_bed_gz_filename(ref_annot_bed_filename):

    ref_annot_bed_gz_filename = ref_annot_bed_filename + ".sorted.bed.gz"
    ref_annot_bed_gz_tabix_filename = ref_annot_bed_gz_filename + ".tbi"

    if not os.path.exists(ref_annot_bed_gz_filename) and not os.path.exists(ref_annot_bed_gz_tabix_filename):
        cmd = f"sort -k1,1 -k2,2g -k3,3g {ref_annot_bed_filename} | bgzip -c > {ref_annot_bed_gz_filename}"
        execute_cmd(cmd)

        execute_cmd(f"tabix {ref_annot_bed_gz_filename}")

    return ref_annot_bed_gz_filename
                        
                        
                        

                        
####
def ensure_expr_matrix_pickle(expr_matrix_filename):

    expr_matrix_pickle = expr_matrix_filename + ".pickle"

    if not os.path.exists(expr_matrix_pickle):
        cmd = " ".join([os.path.join(vif_expr_viewer_basedir, "util/index_expression_matrix.py"),
                                     expr_matrix_filename])
        execute_cmd(cmd)

    return expr_matrix_pickle





if __name__=='__main__':
    main()
