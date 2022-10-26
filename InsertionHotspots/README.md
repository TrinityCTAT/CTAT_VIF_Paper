
# Define insertion hotspots


```    
./scripts/define_insertion_hotspots.py --insertions_tsv data/TCGA.filtered_insertions.tsv --output hotspots

./scripts/plot_hotspots.Rscript --hotspots_tsv CESC.hotspots_defined.tsv
```

Annotate 10 neighboring genes to each side of the insertion 

./scripts/annotate_neighboring_genes.hotspots.py --hotspots hotspots --ref_gene_spans data/ref_annot.gtf.gene_spans.hg38  --num_genes_include 10 --output hotspots.w_neighbors.tsv

    
# Analyze insertion hotspots for CNV and expr outlier enrichments
    
 ./scripts/hotspot_to_expr_n_cnv_analysis.py --ref_gene_spans data/ref_annot.gtf.gene_spans.hg38  --ref_annot_bed data/hg38v22.ref_annot.gtf.bed --expr_matrix data/TCGA_resource_files/TCGA-*fpkm*tsv.gz --cnv_tsv data/TCGA_resource_files/TCGA-*cnv.tsv.gz --hotspots hotspots.w_neighbors.tsv --min_hotspot_samples 2


# Misc example cmds

    ./scripts/hotspot_to_expr_n_cnv_analysis.py --ref_gene_spans data/ref_annot.gtf.gene_spans.hg38  --ref_annot_bed data/hg38v22.ref_annot.gtf.bed --expr_matrix data/TCGA_resource_files/TCGA-*fpkm*tsv.gz --cnv_tsv data/TCGA_resource_files/TCGA-*cnv.tsv.gz --hotspots hotspots.sm_test

    
example vif viewer command:
    ~/GITHUB/CTAT_VIF/VIF_insertion_expression_viewer/vif_expr_viewer.py --ref_annot_bed data/hg38v22.ref_annot.gtf.bed.sorted.bed.gz --ref_gene_spans_bed data/ref_annot.gtf.gene_spans.hg38.bed.gz --expr_matrix_bdbs data/TCGA_resource_files/TCGA-CESC.htseq_fpkm-uq.genesym.tsv.gz.bdb data/TCGA_resource_files/TCGA-LIHC.htseq_fpkm-uq.genesym.tsv.gz.bdb --region chr16:88663298-88716337 --vif_insertions_tsv_tabix_gz hotspots.sm_test.bedlike.tsv.gz --output_filename chr16:88696794.expr_insertions_gview.pdf

    
/Users/bhaas/GITHUB/CTAT_VIF/VIF_insertion_expression_viewer/util/plot_CNV_and_EXPR_outliers_for_cohort_hotspot.heatmap.Rscript --cnv_regions vif.chr16:88663298-88716337.cnv_info.tsv --expr_info vif.chr16:88663298-88716337.expr.tsv --gene_spans vif.chr16:88663298-88716337.gene_spans.tsv --viral_insertions vif.chr16:88663298-88716337.insertions.tsv --output_pdf vif.chr16:88663298-88716337.both_cnv_n_expr.plot.pdf


    

     ~/GITHUB/CTAT_VIF/VIF_insertion_expression_viewer/vif_expr_viewer.py --ref_annot_bed data/hg38v22.ref_annot.gtf.bed.sorted.bed.gz --ref_gene_spans_bed data/ref_annot.gtf.gene_spans.hg38.bed.gz --expr_matrix_bdbs data/TCGA_resource_files/TCGA-CESC.htseq_fpkm-uq.genesym.tsv.gz.bdb data/TCGA_resource_files/TCGA-HNSC.htseq_fpkm-uq.genesyms.tsv.gz.bdb data/TCGA_resource_files/TCGA-LIHC.htseq_fpkm-uq.genesym.tsv.gz.bdb --region chr16:88663298-88716337 --vif_insertions_tsv_tabix_gz hotspots.sm_test.bedlike.tsv.gz --output_filename chr16:88696794.expr_insertions_gview.pdf --cnv_regions_tsv_tabix_gz cnv_info.bedlike.tsv.gz


    
