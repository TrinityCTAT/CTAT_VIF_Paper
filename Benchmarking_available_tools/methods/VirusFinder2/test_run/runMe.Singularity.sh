#!/bin/bash

set -ex


if [ ! -e virusfinder2.simg ]; then
    gsutil cp gs://trinityctatvirusinsertionfinder/SINGULARITY/other_pipelines/virusfinder2.simg .
fi

if [ ! -e subset_insertion_seqs_QSadjust_50ins_7sam_1.fq ]; then

    gsutil cp gs://trinityctatvirusinsertionfinder/OTHER_PIPELINES/VirusFinder2/subset_insertion_seqs_QSadjust_50ins_7sam_1.fq .
    gsutil cp gs://trinityctatvirusinsertionfinder/OTHER_PIPELINES/VirusFinder2/subset_insertion_seqs_QSadjust_50ins_7sam_2.fq .
    
fi


if [ ! -e VIF_virus_reference.tar.gz ]; then
    gsutil cp gs://trinityctatvirusinsertionfinder/OTHER_PIPELINES/VirusFinder2/VIF_virus_reference.tar.gz .
fi

if [ ! -e human_reference.tar.gz ]; then
    gsutil cp gs://trinityctatvirusinsertionfinder/OTHER_PIPELINES/VirusFinder2/human_reference.tar.gz .
fi


if [ ! -d human_reference ]; then
    tar xvf  human_reference.tar.gz
fi

if [ ! -d virus_reference ]; then
    tar xvf VIF_virus_reference.tar.gz
fi

singularity exec -e -W `pwd` -B `pwd`/../ virusfinder2.simg /usr/local/src/VF2Verse_runner.py  --human_resource_dir human_reference --virus_resource_dir virus_reference --target_virus HPV16 --left_fq subset_insertion_seqs_QSadjust_50ins_7sam_1.fq --right_fq subset_insertion_seqs_QSadjust_50ins_7sam_2.fq --sample_id testme


