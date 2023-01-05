#!/bin/bash

if [ ! -d "data_repo" ]; then
    gsutil cp gs://trinityctatvirusinsertionfinder/OTHER_PIPELINES/FastViFi/FastViFi_grch38_data_repo.tar .
    tar xvf FastViFi_grch38_data_repo.tar
fi

if [ ! -d "viral_data" ]; then
    gsutil cp gs://trinityctatvirusinsertionfinder/OTHER_PIPELINES/FastViFi/FastViFi_viral_data_repo.tar .
    tar xvf FastViFi_viral_data_repo.tar
fi

if [ ! -d "kraken_datasets" ]; then
    gsutil cp gs://trinityctatvirusinsertionfinder/OTHER_PIPELINES/FastViFi/kraken_datasets.tar .
    tar xvf kraken_datasets.tar
fi


if [ ! -e "test_reads_1.fq" ]; then
    gsutil cp gs://trinityctatvirusinsertionfinder/OTHER_PIPELINES/FastViFi/test_reads_1.fq .
    gsutil cp gs://trinityctatvirusinsertionfinder/OTHER_PIPELINES/FastViFi/test_reads_2.fq .
fi

if [ ! -e "../WDL/cromwell-58.jar" ]; then
  wget https://github.com/broadinstitute/cromwell/releases/download/58/cromwell-58.jar -O ../WDL/cromwell-58.jar
fi

java -jar ../WDL/cromwell-58.jar \
	run ../WDL/FastViFi.wdl \
	-i  inputs.json
