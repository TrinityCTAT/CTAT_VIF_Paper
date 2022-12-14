#!/bin/bash

set -ex

if [ ! -e hv16.reads_R1.fq ]; then
    gsutil cp gs://trinityctatvirusinsertionfinder/OTHER_PIPELINES/NF_VIF/hv16.reads_R1.fq .
    gsutil cp gs://trinityctatvirusinsertionfinder/OTHER_PIPELINES/NF_VIF/hv16.reads_R2.fq .
fi


if [ ! -e HPV_BWT2_index ]; then
    gsutil cp gs://trinityctatvirusinsertionfinder/OTHER_PIPELINES/NF_VIF/HPV_BWT2_index.tar .
    tar xvf HPV_BWT2_index.tar
fi

if [ ! -e human_reference ]; then
    gsutil cp gs://trinityctatvirusinsertionfinder/OTHER_PIPELINES/NF_VIF/human_reference.tar .
    tar xvf human_reference.tar
fi

if [ ! -e GRCh38.genome.2bit ]; then
    gsutil cp gs://trinityctatvirusinsertionfinder/OTHER_PIPELINES/NF_VIF/GRCh38.genome.2bit .
fi

if [ ! -e ref_genome.fa ]; then
    gsutil cp gs://trinityctatvirusinsertionfinder/OTHER_PIPELINES/NF_VIF/ref_genome.fa .
fi


if [ ! -e "../WDL/cromwell-58.jar" ]; then
  wget https://github.com/broadinstitute/cromwell/releases/download/58/cromwell-58.jar -O ../WDL/cromwell-58.jar
fi

#java -jar ../WDL/cromwell-58.jar \
#	run ../WDL/ViFi.wdl \
#	-i  inputs.json

docker run --rm -v `pwd`:/data trinityctat/nf_vif:devel \
       nextflow /usr/local/src/nf-VIF/main.nf \
       --reads "/data/*_R{1,2}.fq" \
       --genome 'hg38' \
       --bwt2_index /data/human_reference --fasta /data/ref_genome.fa \
       --blatdb /data/GRCh38.genome.2bit  \
       --fasta_hpv /data/HPVs_db.fasta \
       --bwt2_index_hpv  /data/HPV_BWT2_index \
       --outdir /data/hpv16 \
       --split_report  \
       --skip_trimming \
       --skip_fastqc \
       --skip_multiqc \
       --nb_geno 1

