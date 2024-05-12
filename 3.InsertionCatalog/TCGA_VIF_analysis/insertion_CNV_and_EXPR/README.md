

Simplify insertions, make non-redundant:

- first restrict to highest evidence insertion per virus breakpoint group.
- then each sample and within 1Mb region, select the single highest evidence insertion.

```
    ./simplify_best_insertion_per_sample_region.py  --insertions_tsv ../filtered_insertions.tsv --output filtered_insertions.nr.tsv
```


Annotate with neighboring genes and insertion genes, 10 on each side:

```
    ./scripts/annotate_neighboring_genes.insertions.py  --insertions filtered_insertions.nr.tsv --ref_gene_spans ../../InsertionHotspots/data/ref_annot.gtf.gene_spans.hg38 --output filtered_insertions.nr.tsv.w_10neighbors --num_genes_include 10
```


# get expr vals

```
    ./scripts/extract_insertion_expr_info.py  --insertions_tsv filtered_insertions.nr.tsv.w_10neighbors --expr_matrix_bdbs ../../InsertionHotspots/data/TCGA_resource_files/*bdb --output_filename filtered_insertions.nr.tsv.w_10neighbors.expr_vals

```


# get max CNV within region of breakpoint:

```
    ./scripts/extract_insertion_CNV_info.py  --insertions_tsv filtered_insertions.nr.tsv --region_size 10000 --output_filename filtered_insertions.nr.tsv.CNV_within_10kb --cnv_regions_tsv_tabix_gz ../../InsertionHotspots/cnv_info.bedlike.tsv.gz

```    
