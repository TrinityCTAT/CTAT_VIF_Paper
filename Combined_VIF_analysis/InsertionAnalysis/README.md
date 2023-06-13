
# Merge_Virus_Insertions_and_Content.Rmd
- create: all_insertions_and_virus_content_merged.tsv


# Shared_Insertion_Sites.Rmd
- inputs: all_insertions_and_virus_content_merged.tsv
-outputs: insertion_sites_shared_by_diff_participants.tsv
          all_insertions_and_virus_content_merged.w_shared_site_info.tsv (includes 'suspicious' column - includes those w/o matching virus content and those with shared <1% support from dominant sites). Also, defines repr_virus_brkend for representative virus breakend insertions.

# Insertion_threshold_analysis.VirusContentDisparity.Rmd
- examines viral insertions that lack corresponding virus content, and how min read support impacts numbers of insertions and samples found with insertions.
- inputs: all_insertions_and_virus_content_merged.w_shared_site_info.tsv


    
# Insertion_threshold_analysis.by_orthogonal_brkpt_support.Rmd    
- compares RNA-seq based insertions to WXS, WGS, and HYB
- explores minimal read support for supported rna-seq insertions.
- writes: rna_insertions_within_10kb_other_seqtype_insertions.tsv
          
# InsertionCounting.OrthogonalValidations.Rmd
- a follow-up to the above to further examine ortho-validated sites.
- inputs: all_insertions_and_virus_content_merged.w_shared_site_info.tsv  
- outputs: rnaseq_insertions.orthovalidations.tsv

# Filter_Insertions.Rmd
- remove potential viral insertion contaminants
- input: all_insertions_and_virus_content_merged.w_shared_site_info.tsv
- output: all_insertions_and_virus_content_merged.FILTERED.tsv


# InsertionCounting.Rmd
- examine virus insertion counts, w/ or w/o restricting to breakends only.


    
    
##### old stuff:

# Insertion_threshold_analysis.by_orthogonal_brkpt_support.Rmd
- compares RNA-seq based insertions to WXS, WGS, and HYB
- explores minimal read support for supported rna-seq insertions.
- writes: rna_insertions_within_10kb_other_seqtype_insertions.tsv  ## - note this is read by others above!
        
# Insertion_threshold_analysis.VirusContentDisparity.Rmd
- finds insertions where there is no robust support for matching virus content
- writes: rna_insertions_and_virus_content_merged.tsv and rna_insertions_without_matching_virus_content.tsv

        
# Insertions.SupsiciousSharedSites.Rmd
- examine potential source of virus insertion derived from bleed-thru sequeencing.
- writes: bleedthru_insertion_candidates.tsv

