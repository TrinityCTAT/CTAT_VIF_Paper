#!/bin/bash

set -ex

if [ ! -e HPV16.noCHR_left.fq ]; then
    gsutil cp gs://trinityctatvirusinsertionfinder/OTHER_PIPELINES/BatVI/HPV16.noCHR_left.fq .
    gsutil cp gs://trinityctatvirusinsertionfinder/OTHER_PIPELINES/BatVI/HPV16.noCHR_right.fq .
fi


if [ ! -e BLAST_index_noCHR.tar.gz ]; then
    gsutil cp gs://trinityctatvirusinsertionfinder/OTHER_PIPELINES/BatVI/BLAST_index_noCHR.tar.gz .
fi

if [ ! -e Human_BatIndex.tar.gz ]; then
    gsutil cp gs://trinityctatvirusinsertionfinder/OTHER_PIPELINES/BatVI/Human_BatIndex.tar.gz .
fi

if [ ! -e BWA_index_noCHR.tar.gz ]; then
    gsutil cp gs://trinityctatvirusinsertionfinder/OTHER_PIPELINES/BatVI/BWA_index_noCHR.tar.gz .
fi

if [ ! -e HPV_BLAST_index.tar.gz ]; then
    gsutil cp gs://trinityctatvirusinsertionfinder/OTHER_PIPELINES/BatVI/HPV_BLAST_index.tar.gz .
fi


if [ ! -e Virus_BatIndex.tar.gz ]; then
    gsutil cp gs://trinityctatvirusinsertionfinder/OTHER_PIPELINES/BatVI/Virus_BatIndex.tar.gz .
fi

if [ ! -e HPV_BWA_index.tar.gz ]; then
    gsutil cp gs://trinityctatvirusinsertionfinder/OTHER_PIPELINES/BatVI/HPV_BWA_index.tar.gz .
fi


if [ ! -e "../WDL/cromwell-58.jar" ]; then
  wget https://github.com/broadinstitute/cromwell/releases/download/58/cromwell-58.jar -O ../WDL/cromwell-58.jar
fi

java -jar ../WDL/cromwell-58.jar \
	run ../WDL/BatVI.wdl \
	-i  inputs.json

