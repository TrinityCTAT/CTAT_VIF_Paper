#!/bin/bash

singularity build  virusfinder2.simg  docker://trinityctat/virusfinder2:devel 


singularity exec -e -C virusfinder2.simg /usr/local/src/VF2Verse_runner.py --help

# copy to:
# gs://trinityctatvirusinsertionfinder/SINGULARITY/other_pipelines/virusfinder2.simg
