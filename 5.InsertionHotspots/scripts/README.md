# define hotspots
define_insertion_hotspots.py --insertions_tsv CESC.refined_results.agg.tsv.gz --output CESC.hotspots_defined.tsv

# plot hotspot manhattan style
plot_hotspots.Rscript --hotspots_tsv CESC.hotspots_defined.tsv


# add neighboring genes
annotate_neighboring_genes.hotspots.py --hotspots CESC.hotspots_defined.tsv --ref_gene_spans ref_annot.gtf.gene_spans.hg38 --output CESC.hotspots_defined.w_neighbors.tsv


