# Order of operations for Combining data and analyses.

## Combine separately processed cohort results:

    
Combine_Processed_Samples_Data.Rmd :  generate 'all_samples_processed.tsv'

Combine_Virus_Content_Data.Rmd : generate 'virus_content_unfiltered.tsv'

Combine_Virus_Insertions_Data.Rmd : generate 'virus_insertions_unfiltered.tsv'


## Merge into a single insertion and content file

InsertionAnalysis/Merge_Virus_Insertions_and_Content.Rmd : generates 'all_insertions_and_virus_content_merged.tsv', merging 'virus_insertions_unfiltered.tsv' and 'virus_content_unfiltered.tsv', and incorporating the 'matching_virus_content' boolean column.

### Examine shared insertion sites (potential contaminants)

InsertionAnalysis/Shared_Insertion_Sites.Rmd :

    Adds a 'suspicious' column, if the insertion lacks matching virus content or is a shared site that's less than 1% the dominant read support of another shared site.
    - input: all_insertions_and_virus_content_merged.tsv
    - generates 'all_insertions_and_virus_content_merged.w_shared_site_info.tsv', 

### Filter out suspicious sites:

InsertionAnalysis/Filter_Insertions.Rmd :

    Excludes the suspicious entries:
    - input: all_insertions_and_virus_content_merged.w_shared_site_info.tsv
    - output: all_insertions_and_virus_content_merged.FILTERED.tsv
    

## Compute Insertion Hotspots:

    HotspotsRevisited/
    - run: ./run_hotspots.sh
       - input: ../InsertionAnalysis/all_insertions_and_virus_content_merged.FILTERED.tsv
       - output: hotspots.win_1e5.tsv.w_20_neighbors.regrouped_by_insert_gene
                 plot.hotspots.win_1e5.tsv.w_20_neighbors.regrouped_by_insert_gene.hotspot_virus_sample_counts.tsv


## Decorate Insertions with hotspot, fusions, and expression results.

    DecorateInsertions/DecorateInsertions.Rmd

    -inputs: ../HotspotsRevisited/hotspots.win_1e5.tsv.w_20_neighbors.regrouped_by_insert_gene
             ../Insertion_STARF_Fusion_Mapping/FI_fusions_at_insertions.tsv
             ../Insertion_Spliced_Human/spliced_insertions.tsv
             ../Insertion_and_CNVs/TCGA_insertions_and_CNV_within10kb
             ../Insertion_and_EXPR_effects/virus_insertions_unfiltered.10kb.expression_region_analysis.tsv


    -incorporates columns: fusions_at_insertion, spliced_gene, copy_number, and expr_quantile
    
    -output: all_insertions_and_virus_content_merged.FILTERED.Decorated.tsv
    
## Examine the hotspots:

### Evaluate Size Distributions for Random Hotspot Sizes

    HotspotsRevisited/AnalyzeHotspotSizes.Rmd
    
### Examine the top hotspots, virus and tissue type representation according to cohort and seqtypes

    HotspotsRevisited/AnalyzeTopHotspots.Rmd



    

    
# Misc analyses

## Benchmarking ctat-VIF with simulated data

    Fig 1:
    
    characterizing_VIF_accuracy_wSIM_data/uniform_coverage_virus_detection/PE/sim50_PE/compare_prelim_to_refined.Rmd : Figures for sensitivity barchars and  precision / recall curve 
    characterizing_VIF_accuracy_wSIM_data/variable_coverage_virus_detection/QuantComparisonPrelimVsRefined.Rmd : Figure for evidence read quantification agreement 




    
## Benchmarking against other tools

    Figure 2:
    
    Benchmarking_available_tools/__bmark_bhaas_sim50_WITHREPS/examine_bmark_accuracy.by_replicate.Rmd
    
    
## Comparison to NATGEN2015 study results

    Supp Figure 1:

    NatGenetics2015_reanalysis/NatGenetics2015_VIF/compare_VIF_vs_published_insertions/compare_VIF_to_published_NatGenet2015.Rmd
    NatGenetics2015_reanalysis/NatGenetics2015_VIF/compare_VIF_vs_published_insertions/compare_VIF_to_published_NatGenet2015.SangerValidated.Rmd


# Application to full collection of samples:

    Supp Figure 2:
    
    Combined_VIF_analysis/SampleCounting.Rmd


## Virus content

    Figure 3
    Supp Figure 3

    Combined_VIF_analysis/VirusContentAnalysis/VirusContentAnalysis.Rmd


## Insertion content

    Supp. Figure 4
    Combined_VIF_analysis/InsertionAnalysis/Insertion_threshold_analysis.VirusContentDisparity.Rmd

    Combined_VIF_analysis/InsertionAnalysis/Shared_Insertion_Sites.Rmd

    Combined_VIF_analysis/InsertionAnalysis/InsertionCounting.OrthogonalValidations.Rmd
    Combined_VIF_analysis/InsertionAnalysis/Insertion_threshold_analysis.by_orthogonal_brkpt_support.Rmd

    Supp. Figure 5
    Combined_VIF_analysis/InsertionAnalysis/InsertionCounting.Rmd


## Functional impacts:

    Combined_VIF_analysis/Insertion_and_CNVs/Insertions_and_CNVs.Rmd

    Combined_VIF_analysis/Insertion_and_EXPR_effects/InsertionExpressionEffects.Rmd
    
    Combined_VIF_analysis/Insertion_STARF_Fusion_Mapping/InsertionFusionMappingCurrentFusionInspector.Rmd
