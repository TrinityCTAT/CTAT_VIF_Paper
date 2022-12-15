#!/bin/bash

set -ex

ctat-vif --genome_lib_dir /data \
         --viral_fasta /data/virus_db.fasta \
         --star_index_plus_virus /data/hg_plus_viraldb.fasta.star.idx \
         --cpu 4 \
         --left /data/vif_reads_1.fastq.gz \
         --right /data/vif_reads_2.fastq.gz \
         --min_reads 5 \
         --max_hits 50 \
         --sample_id test \
         -O /data/vif_outdir


