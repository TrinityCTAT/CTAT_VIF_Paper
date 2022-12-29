all: ctat_vif vifi VIRUSBreakend BatVI VF2Verse NF_VIF



ctat_vif:
	../bmark_scripts/aggregate_ctat_vif.py data/ctat_vif/*vif.refined.tsv > data/ctat_vif.agg.tsv
	../bmark_scripts/prep_ctat_vif_for_bmarking.py data/ctat_vif.agg.tsv  > data/ctat_vif.agg.tsv.prepped
	../bmark_scripts/map_to_truth_insertions.py --truth_insertions ./insertion_truth_set.tsv --predicted_insertions data/ctat_vif.agg.tsv.prepped  > data/ctat_vif.agg.tsv.prepped.truthset_mapped
	../bmark_scripts/analyze_VIF_TP_FP_FNs.Rscript --dat data/ctat_vif.agg.tsv.prepped.truthset_mapped


vifi:
	../bmark_scripts/aggregate_vifi.py data/vifi/*clusters.txt > data/vifi.agg.tsv
	../bmark_scripts/prep_ViFi_for_bmarking.py data/vifi.agg.tsv > data/vifi.agg.tsv.prepped
	../bmark_scripts/map_to_truth_insertions.py --truth_insertions ./insertion_truth_set.tsv --predicted_insertions data/vifi.agg.tsv.prepped   > data/vifi.agg.tsv.prepped.truthset_mapped      
	../bmark_scripts/analyze_VIF_TP_FP_FNs.Rscript --dat data/vifi.agg.tsv.prepped.truthset_mapped



VIRUSBreakend:
	../bmark_scripts/aggregate_VIRUSBreakend.py data/virusbreakend/*.vcf > data/virusbreakend.agg.tsv
	../bmark_scripts/prep_virusbreakend_for_benchmarking.py data/virusbreakend.agg.tsv > data/virusbreakend.agg.tsv.prepped
	../bmark_scripts/map_to_truth_insertions.py --truth_insertions ./insertion_truth_set.tsv --predicted_insertions data/virusbreakend.agg.tsv.prepped   > data/virusbreakend.agg.tsv.prepped.truthset_mapped      
	../bmark_scripts/analyze_VIF_TP_FP_FNs.Rscript --dat data/virusbreakend.agg.tsv.prepped.truthset_mapped



BatVI:
	../bmark_scripts/aggregate_BatVI.py data/batvi/*final_hits.txt > data/batvi.agg.tsv
	../bmark_scripts/prep_batvi_for_bmarking.py data/batvi.agg.tsv > data/batvi.agg.tsv.prepped
	../bmark_scripts/map_to_truth_insertions.py --truth_insertions ./insertion_truth_set.tsv --predicted_insertions data/batvi.agg.tsv.prepped   > data/batvi.agg.tsv.prepped.truthset_mapped      
	../bmark_scripts/analyze_VIF_TP_FP_FNs.Rscript --dat data/batvi.agg.tsv.prepped.truthset_mapped


VF2Verse:
	../bmark_scripts/aggregate_VF2Verse.py data/VF2Verse/*.txt > data/VF2Verse.agg.tsv
	../bmark_scripts/prep_VF2Verse_for_bmarking.py data/VF2Verse.agg.tsv > data/VF2Verse.agg.tsv.prepped
	../bmark_scripts/map_to_truth_insertions.py --truth_insertions ./insertion_truth_set.tsv --predicted_insertions data/VF2Verse.agg.tsv.prepped   > data/VF2Verse.agg.tsv.prepped.truthset_mapped      
	../bmark_scripts/analyze_VIF_TP_FP_FNs.Rscript --dat data/VF2Verse.agg.tsv.prepped.truthset_mapped


NF_VIF:
	../bmark_scripts/aggregate_NF_VIF.py data/nf_vif/*nf-vif > data/nf_vif.agg.tsv
	../bmark_scripts/prep_NF_VIF_for_bmarking.py data/nf_vif.agg.tsv > data/nf_vif.agg.tsv.prepped
	../bmark_scripts/map_to_truth_insertions.py --truth_insertions ./insertion_truth_set.tsv --predicted_insertions  data/nf_vif.agg.tsv.prepped   >  data/nf_vif.agg.tsv.prepped.truthset_mapped
	../bmark_scripts/analyze_VIF_TP_FP_FNs.Rscript --dat data/nf_vif.agg.tsv.prepped.truthset_mapped


