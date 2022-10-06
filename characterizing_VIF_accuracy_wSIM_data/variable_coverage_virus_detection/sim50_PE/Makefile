

all:
# Assign each predicted insertion to its closest truth insertion site
	../../bin/examine_sim_insertions.py  --truth_insertions data/insertion_truth_set.tsv.gz  --predicted_insertions data/prelim_results.agg.tsv.gz > prelim_insertion_mappings.tsv
	../../bin/examine_sim_insertions.py  --truth_insertions data/insertion_truth_set.tsv.gz  --predicted_insertions data/refined_results.agg.tsv.gz > refined_insertion_mappings.tsv 

# Perform accuracy analysi
	../../bin/analyze_VIF_TP_FP_FNs.Rscript --dat prelim_insertion_mappings.tsv --score_vir_brkend_grp_only
	../../bin/analyze_VIF_TP_FP_FNs.Rscript --dat refined_insertion_mappings.tsv --score_vir_brkend_grp_only


clean:
	rm -f ./prelim_insertion_mappings.tsv*
	rm -f ./refined_insertion_mappings.tsv*

