test: prep_inputs
	../vif_expr_viewer.py --ref_annot_bed ref_annot.bed.gz --ref_gene_spans_bed_tabix_gz ref_annot.gtf.gene_spans.bed.gz --expr_matrix_pickle chr8^128500000.expr.matrix.pickle --region chr8:127000000-129000000 --vif_insertions_tsv_tabix_gz hotspot.refined.sample_data.tsv.bedlike.tsv.gz --gistic_matrix_pickle gistic.matrix.pickle --cnv_regions_tsv_tabix_gz tcga_gdc.CNV_info.tsv.bedlike.tsv.gz  --output_filename test_example.pdf



prep_inputs:
#   prep gene spans bed
	../util/gtf_gene_spans_to_bed.py data/ref_annot.gtf.gene_spans | sort -k 1,1 -k2,2g -k3,3g | bgzip -c >   ref_annot.gtf.gene_spans.bed.gz && tabix ref_annot.gtf.gene_spans.bed.gz
#   prep ref annot bed
	cat data/ref_annot.bed | sort -k 1,1 -k2,2g -k3,3g | bgzip -c > ref_annot.bed.gz && tabix ref_annot.bed.gz
#   prep gene expr matrix
	ln -sf data/chr8^128500000.expr.matrix && ../util/index_expression_matrix.py chr8^128500000.expr.matrix
#   prep hotspot data summary
	ln -sf data/hotspot.refined.sample_data.tsv && ../util/hotspots_to_bedlike_tsv.py hotspot.refined.sample_data.tsv && bgzip -f hotspot.refined.sample_data.tsv.bedlike.tsv && tabix -p bed hotspot.refined.sample_data.tsv.bedlike.tsv.gz
# prep gistic
	ln -sf data/gistic.matrix && ../util/index_expression_matrix.py gistic.matrix
# prep cnv regions info
	ln -sf data/tcga_gdc.CNV_info.tsv && ../util/cnv_to_bedlike_tsv.py tcga_gdc.CNV_info.tsv && bgzip -f tcga_gdc.CNV_info.tsv.bedlike.tsv && tabix -p bed tcga_gdc.CNV_info.tsv.bedlike.tsv.gz
	touch prep_inputs

clean:
	rm -f ./ref_annot.*
	rm -f ./hotspot.*
	rm -f ./chr8*
	rm -f ./test_example.pdf
	rm -f ./gistic.*
	rm -f ./tcga_gdc.CNV_info.*
	rm -f ./prep_inputs



just_cnv_plot:
	../util/plot_CNV_insertions_for_cohort.Rscript --cnv_regions data/cnv.chr8:127000000-129000000.cnv_info.tsv --gene_spans data/cnv.chr8:127000000-129000000.gene_spans.tsv --viral_insertions data/cnv.chr8:127000000-129000000.insertions.tsv --output_pdf test.cnv.pdf

just_expr_plot:
	../util/plot_EXPR_outliers_for_cohort_hotspot.heatmap.Rscript  --expr_info data/cnv.chr8\:127000000-129000000.expr.tsv --gene_spans data/cnv.chr8\:127000000-129000000.gene_spans.tsv --viral_insertions data/cnv.chr8\:127000000-129000000.insertions.tsv --output_pdf test.expr.pdf


just_combined_cnv_expr_insertion_plot:
	../util/plot_CNV_and_EXPR_outliers_for_cohort_hotspot.heatmap.Rscript  --cnv_regions data/cnv.chr8:127000000-129000000.cnv_info.tsv --expr_info data/cnv.chr8\:127000000-129000000.expr.tsv --gene_spans data/cnv.chr8\:127000000-129000000.gene_spans.tsv --viral_insertions data/cnv.chr8\:127000000-129000000.insertions.tsv --output_pdf test.expr_n_cnv.pdf


just_expr_wilcoxon_plot:
	../util/plot_EXPR_outliers_for_cohort_hotspot.gene_Wilcoxons.Rscript  --expr_info data/cnv.chr8\:127000000-129000000.expr.tsv --gene_spans data/cnv.chr8\:127000000-129000000.gene_spans.tsv --viral_insertions data/cnv.chr8\:127000000-129000000.insertions.tsv --output_pdf test.expr.wilcoxon.pdf

just_cnv_wilcoxon_plot:
	../util/plot_CNV_for_cohort_hotspot.gene_Wilcoxons.Rscript --cnv_regions data/cnv.chr8:127000000-129000000.cnv_info.tsv --gene_spans data/cnv.chr8:127000000-129000000.gene_spans.tsv --viral_insertions data/cnv.chr8:127000000-129000000.insertions.tsv --output_pdf test.cnv.wilcoxon.pdf

