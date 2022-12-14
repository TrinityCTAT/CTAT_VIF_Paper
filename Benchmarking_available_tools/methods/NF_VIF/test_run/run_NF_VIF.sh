#!/bin/bash

set -ex

cd /data

nextflow /usr/local/src/nf-VIF/main.nf \
         --reads "/data/*_R{1,2}.fq" \
         --genome hg38 \
         --bwt2_index /data/human_reference \
         --fasta /data/ref_genome.fa \
         --blatdb /data/GRCh38.genome.2bit \
         --fasta_hpv /data/HPVs_db.fasta \
         --bwt2_index_hpv  /data/HPV_BWT2_index \
         --outdir /data/hpv16 \
         --split_report --skip_trimming --skip_fastqc --skip_multiqc --nb_geno 1 

