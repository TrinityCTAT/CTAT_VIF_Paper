#!/bin/bash

set -ex


WINSIZE=1e5
WINSIZE_INT=100000
MIN_HOTSPOT_SAMPLES=3
MIN_EV_READS=1


scripts/define_insertion_hotspots.py \
    --insertions_tsv ../InsertionAnalysis/all_insertions_and_virus_content_merged.FILTERED.tsv \
    --window_size $WINSIZE_INT \
    --min_ev_reads $MIN_EV_READS \
    --output hotspots.win_$WINSIZE.tsv


scripts/annotate_neighboring_genes.hotspots.py --hotspots hotspots.win_$WINSIZE.tsv  --ref_gene_spans ../../InsertionHotspots/data/ref_annot.gtf.gene_spans.hg38 --num_genes_include 20 --output hotspots.win_$WINSIZE.tsv.w_20_neighbors --no_gene_decoration


scripts/regroup_hotspots_by_gene.Rscript hotspots.win_$WINSIZE.tsv.w_20_neighbors


scripts/plot_hotspots.v2.Rscript --hotspots_tsv hotspots.win_$WINSIZE.tsv.w_20_neighbors.regrouped_by_insert_gene --min_hotspot_samples $MIN_HOTSPOT_SAMPLES --label_hotspot_min 11


