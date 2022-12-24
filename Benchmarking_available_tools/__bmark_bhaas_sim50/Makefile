all: ctat_vif vifi VIRUSBreakend BatVI



ctat_vif:
	../bmark_scripts/aggregate_ctat_vif.py data/ctat_vif/*vif.refined.tsv > data/ctat_vif.agg.tsv

vifi:
	../bmark_scripts/aggregate_vifi.py data/vifi/*clusters.txt > data/vifi.agg.tsv


VIRUSBreakend:
	../bmark_scripts/aggregate_VIRUSBreakend.py data/virusbreakend/*.vcf > data/virusbreakend.agg.tsv


BatVI:
	../bmark_scripts/aggregate_BatVI.py data/batvi/*final_hits.txt > data/batvi.agg.tsv



