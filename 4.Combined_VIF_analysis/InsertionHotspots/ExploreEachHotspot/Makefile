

make_gview_plots:
	../scripts/hotspot_to_expr_n_cnv_analysis.py --ref_gene_spans ../../../data/ref_annot.gtf.gene_spans.hg38.gz --ref_annot_bed ../../../data/hg38v22.ref_annot.gtf.bed.sorted.bed.gz --expr_matrix TCGA_resource_files/TCGA-*fpkm*tsv.gz --cnv_tsv TCGA_resource_files/TCGA-*cnv.tsv.gz --hotspots ../hotspots.win_1e5.tsv.w_20_neighbors.regrouped_by_insert_gene.gz --min_hotspot_samples 11 --expr_bed_dir  TCGA_resource_files


prep_make_all_misc_plots:
	~/GITHUB/CTAT_VIF/VIF_insertion_expression_viewer/hotspot_to_misc_plots.py --hotspots_bed_gz ../hotspots.win_1e5.tsv.w_20_neighbors.regrouped_by_insert_gene.gz.bedlike.tsv.gz --expr_bed_dir TCGA_resource_files/ --min_hotspot_samples 11 --cnv_bed_gz cnv_info.bedlike.tsv.gz --ref_annot_bed_gz ../../../data/hg38v22.ref_annot.gtf.bed.sorted.bed.gz --ref_gene_spans_bed_gz ../../../data/ref_annot.gtf.gene_spans.hg38.gz.bed.gz

TCGA_resource_files:
	wget --no-parent -nH -r --cut-dirs=2  https://data.broadinstitute.org/Trinity/VIF_PAPER/TCGA_resource_files/ 
#   as backup, also here: gs://ctat-virusintegrationfinder-archive/VIF_PAPER/TCGA_resource_files.tar.gz

clean:
	rm -rf ./TCGA_resource_files
