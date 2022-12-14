#!/bin/bash

set -ex

if [ ! -e reads_R1.fq ]; then
    gsutil cp gs://trinityctatvirusinsertionfinder/OTHER_PIPELINES/NF_VIF/hpv16.reads_R1.fq reads_R1.fq
    gsutil cp gs://trinityctatvirusinsertionfinder/OTHER_PIPELINES/NF_VIF/hpv16.reads_R2.fq reads_R2.fq
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


if [ ! -e HPVs_db.fasta  ]; then
    gsutil cp gs://trinityctatvirusinsertionfinder/OTHER_PIPELINES/NF_VIF/HPVs_db.fasta .
fi


#if [ ! -e "../WDL/cromwell-58.jar" ]; then
#  wget https://github.com/broadinstitute/cromwell/releases/download/58/cromwell-58.jar -O ../WDL/cromwell-58.jar
#fi

#java -jar ../WDL/cromwell-58.jar \
#	run ../WDL/ViFi.wdl \
#	-i  inputs.json


docker run -e USER='myself' --rm -it -v `pwd`:/data -v /tmp:/tmp  trinityctat/nf_vif:devel bash /data/run_NF_VIF.sh

#docker run -e USER='myself' --rm -it -v `pwd`:/data -v /tmp:/tmp   brownmp/nf_vif:devel bash /data/run_NF_VIF.sh

