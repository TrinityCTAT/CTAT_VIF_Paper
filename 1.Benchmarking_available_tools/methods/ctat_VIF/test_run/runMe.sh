#!/bin/bash

set -ex

if [ ! -e vif_reads_1.fastq.gz ]; then
    gsutil cp gs://trinityctatsampledata/ctat-virusintegrationfinder/vif_reads_1.fastq.gz .
    gsutil cp gs://trinityctatsampledata/ctat-virusintegrationfinder/vif_reads_2.fastq.gz .
fi



if [ ! -e ref_genome.fa ]; then
    gsutil cp gs://ctat_genome_libs/GRCh38_gencode_v22/03-01-2021/ref_genome.fa .
fi

if [ ! -e ref_annot.gtf ]; then
    gsutil cp gs://ctat_genome_libs/GRCh38_gencode_v22/03-01-2021/ref_annot.gtf .
fi

if [ ! -e ref_genome.fa.star.idx ]; then
    gsutil cp gs://ctat_genome_libs/GRCh38_gencode_v22/03-01-2021/ref_genome.fa.star.idx.tar .
    tar xvf ref_genome.fa.star.idx.tar
    rm -f ref_genome.fa.star.idx.tar
fi


if [ ! -e hg_plus_viraldb.fasta.star.idx ]; then
    gsutil cp gs://ctat_genome_libs/GRCh38_gencode_v22/06-08-2022/hg_plus_viraldb.fasta.star.idx.tar .
    tar xvf hg_plus_viraldb.fasta.star.idx.tar
    rm -f hg_plus_viraldb.fasta.star.idx.tar
fi


if [ ! -e virus_db.fasta ]; then
    gsutil cp gs://ctat_genome_libs/GRCh38_gencode_v22/06-08-2022/virus_db.fasta .
fi


# run via docker

docker run --rm -it -v /tmp:/tmp -v `pwd`:/data trinityctat/ctat_vif:latest bash /data/run_ctat-VIF.sh


