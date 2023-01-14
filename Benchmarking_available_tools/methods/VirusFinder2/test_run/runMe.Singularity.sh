#!/bin/bash

set -ex


if [ ! -e virusfinder2.simg ]; then
    gsutil cp gs://trinityctatvirusinsertionfinder/SINGULARITY/other_pipelines/virusfinder2.simg .
fi

if [ ! -e rep1.HPV33.insertion_seqs_1.fq.gz ]; then

    gsutil cp gs://trinityctatvirusinsertionfinder/OTHER_PIPELINES/VirusFinder2/rep1.HPV33.insertion_seqs_1.fq.gz .
    gsutil cp gs://trinityctatvirusinsertionfinder/OTHER_PIPELINES/VirusFinder2/rep1.HPV33.insertion_seqs_2.fq.gz .
    
    #gsutil cp gs://trinityctatvirusinsertionfinder/OTHER_PIPELINES/VirusFinder2/subset_insertion_seqs_QSadjust_50ins_7sam_1.fq .
    #gsutil cp gs://trinityctatvirusinsertionfinder/OTHER_PIPELINES/VirusFinder2/subset_insertion_seqs_QSadjust_50ins_7sam_2.fq .
    
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

if [ ! -d singularity_outdir ]; then
    mkdir singularity_outdir
fi


rootdir=`pwd`

cd singularity_outdir && cp ../*.fq.gz . && singularity exec -e -B $rootdir  ../virusfinder2.simg /usr/local/src/VF2Verse_runner.py  --human_resource_dir `pwd`/../human_reference --virus_resource_dir `pwd`/../virus_reference --target_virus HPV33 --left_fq rep1.HPV33.insertion_seqs_1.fq.gz --right_fq rep1.HPV33.insertion_seqs_2.fq.gz --sample_id testme


