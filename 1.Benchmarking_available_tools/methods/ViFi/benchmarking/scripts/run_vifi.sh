#!/bin/bash


###
# Create the  reference and index 
#cat HG_Resources/GRCh38.genome.fa Virus_DB/all/HPV_plus.fasta > Virus_DB/all/GRCh38_all.fas
#bwa index Virus_DB/all/GRCh38_all.fas



sudo rm -rf cromwell-*


java -jar /home/mbrown/Tools/cromwell-71.jar \
	run /home/mbrown/GitHub/Maxwell_ViFi/ViFi/WDL/ViFi.wdl \
	-i  inputs.json
