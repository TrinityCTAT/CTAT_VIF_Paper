]1;95;0c#!/bin/bash

set -ex

if [ ! -e HPV16.left.fq ]; then
    gsutil cp gs://trinityctatvirusinsertionfinder/OTHER_PIPELINES/VIRUSBreakend/HPV16.left.fq .
    gsutil cp gs://trinityctatvirusinsertionfinder/OTHER_PIPELINES/VIRUSBreakend/HPV16.right.fq .
fi


if [ ! -e virusbreakenddb_20210401.tar.gz ]; then
    gsutil cp gs://trinityctatvirusinsertionfinder/OTHER_PIPELINES/VIRUSBreakend/virusbreakenddb_20210401.tar.gz .
fi

if [ ! -e BWA_index.tar.gz ]; then
    gsutil cp gs://trinityctatvirusinsertionfinder/OTHER_PIPELINES/VIRUSBreakend/BWA_index.tar.gz .
fi


if [ ! -e "../WDL/cromwell-58.jar" ]; then
  wget https://github.com/broadinstitute/cromwell/releases/download/58/cromwell-58.jar -O ../WDL/cromwell-58.jar
fi

java -jar ../WDL/cromwell-58.jar \
	run ../WDL/VIRUSBreakend.wdl \
	-i  inputs.json


