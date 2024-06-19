

make_all_plots:
	../scripts/hotspot_to_expr_n_cnv_analysis.py --ref_gene_spans ../../../data/ref_annot.gtf.gene_spans.hg38.gz --ref_annot_bed ../../../data/hg38v22.ref_annot.gtf.bed.sorted.bed.gz --expr_matrix TCGA_resource_files/TCGA-*fpkm*tsv.gz --cnv_tsv TCGA_resource_files/TCGA-*cnv.tsv.gz --hotspots ../hotspots.win_1e5.tsv.w_20_neighbors.regrouped_by_insert_gene.gz --min_hotspot_samples 11 --expr_bed_dir  TCGA_resource_files


TCGA_resource_files:
	wget --no-parent -nH -r --cut-dirs=2  https://data.broadinstitute.org/Trinity/VIF_PAPER/TCGA_resource_files/ 


clean:
	rm -rf ./TCGA_resource_files
