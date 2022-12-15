all:
	./runMe.sh


cleanout:
	sudo rm -rf ./vif_outdir ./cromwell-executions/ ./cromwell-workflow-logs/


clean: cleanout
	rm -f ./vif_reads_1.fastq.gz ./vif_reads_2.fastq.gz ./ref_genome.fa ./ref_annot.gtf  ./virus_db.fasta
	rm -rf ./ref_genome.fa.star.idx
	rm -rf ./hg_plus_viraldb.fasta.star.idx

