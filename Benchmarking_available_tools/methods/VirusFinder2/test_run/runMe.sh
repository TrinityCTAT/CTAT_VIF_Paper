#!/bin/bash

set -ex

if [ ! -e subset.insertion_seqs_QSadjust_50ins_7sam_left.fq ]; then
    gsutil cp gs://trinityctatvirusinsertionfinder/OTHER_PIPELINES/VIFI/subset.insertion_seqs_QSadjust_50ins_7sam_left.fq .
    gsutil cp gs://trinityctatvirusinsertionfinder/OTHER_PIPELINES/VIFI/subset.insertion_seqs_QSadjust_50ins_7sam_right.fq .
fi


if [ ! -e VIF_virus_reference.tar.gz ]; then
    gsutil cp gs://trinityctatvirusinsertionfinder/OTHER_PIPELINES/VirusFinder2/VIF_virus_reference.tar.gz .
fi

if [ ! -e human_reference.tar.gz ]; then
    gsutil cp gs://trinityctatvirusinsertionfinder/OTHER_PIPELINES/VirusFinder2/human_reference.tar.gz .
fi

if [ ! -e "../WDL/cromwell-58.jar" ]; then
  wget https://github.com/broadinstitute/cromwell/releases/download/58/cromwell-58.jar -O ../WDL/cromwell-58.jar
fi

java -jar ../WDL/cromwell-58.jar \
	run ../WDL/VirusFinder2 \
	-i  inputs.json


