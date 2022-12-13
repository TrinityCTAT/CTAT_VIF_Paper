#!/bin/bash

if [ ! -e subset.insertion_seqs_QSadjust_50ins_7sam_left.fq ]; then
    gsutil cp gs://trinityctatvirusinsertionfinder/OTHER_PIPELINES/VIFI/subset.insertion_seqs_QSadjust_50ins_7sam_left.fq .
    gsutil cp gs://trinityctatvirusinsertionfinder/OTHER_PIPELINES/VIFI/subset.insertion_seqs_QSadjust_50ins_7sam_right.fq .
fi


if [ ! -e Virus_DB.tar.gz ]; then
    gsutil cp gs://trinityctatvirusinsertionfinder/OTHER_PIPELINES/VIFI/Virus_DB.tar.gz .
fi

if [ ! -e Human_DB.tar.gz ]; then
    gsutil cp gs://trinityctatvirusinsertionfinder/OTHER_PIPELINES/VIFI/Human_DB.tar.gz .
fi

