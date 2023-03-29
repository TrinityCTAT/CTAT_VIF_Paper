#!/bin/bash

set -ex


WINSIZE=4e5
WINSIZE_INT=400000
MIN_HOTSPOT_SAMPLES=3


../../InsertionHotspots/scripts/define_insertion_hotspots.py --insertions_tsv all_insertions_and_virus_content_merged.FILTERED.tsv.distilled.1e+07 --window_size $WINSIZE_INT --output hotspots.win_$WINSIZE.tsv


../../InsertionHotspots/scripts/annotate_neighboring_genes.hotspots.py --hotspots hotspots.win_$WINSIZE.tsv  --ref_gene_spans ../../InsertionHotspots/data/ref_annot.gtf.gene_spans.hg38 --num_genes_include 20 --output hotspots.win_$WINSIZE.tsv.w_20_neighbors --no_gene_decoration

../../InsertionHotspots/scripts/plot_hotspots.v2.Rscript --hotspots_tsv hotspots.win_$WINSIZE.tsv.w_20_neighbors --min_hotspot_samples $MIN_HOTSPOT_SAMPLES


