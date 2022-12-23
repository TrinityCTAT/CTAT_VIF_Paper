all: ctat_vif vifi VIRUSBreakend



ctat_vif:
	../bmark_scripts/aggregate_ctat_vif.py data/ctat_vif/*vif.refined.tsv > data/ctat_vif.agg.tsv

vifi:
	../bmark_scripts/aggregate_vifi.py data/vifi/*clusters.txt > data/vifi.agg.tsv


VIRUSBreakend:
	../bmark_scripts/aggregate_VIRUSBreakend.py data/virusbreakend/*.vcf > data/virusbreakend.agg.tsv



