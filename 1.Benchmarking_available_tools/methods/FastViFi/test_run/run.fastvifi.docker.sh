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



docker run --rm -v `pwd`:/data -e VIFI_DIR=/usr/local/src/ViFi -e REFERENCE_REPO=/data/viral_data -e AA_DATA_REPO=/data/data_repo trinityctat/fastvifi:devel python /usr/local/src/FastViFi/run_kraken_vifi_pipeline.py --output-dir /data/outdir --input-file /data/test_reads_1.fq --input-file-2 /data/test_reads_2.fq --level sample-level-validation-intermediate --kraken-path /usr/local/bin/kraken2 --kraken-db-path /data/kraken_datasets --vifi-path /usr/local/src/ViFi/scripts/run_vifi.py --virus hpv --human-chr-list /usr/local/src/FastViFi/test/human_chr_list.txt --skip-bwa-filter --keep-intermediate-files --vifi-human-ref-dir /data/data_repo --vifi-viral-ref-dir /data/viral_data --docker

