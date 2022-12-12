#!/bin/bash

set -e

# if want to remove the previous runs 
sudo rm -rf cromwell-*


java -jar /home/mbrown/Tools/cromwell-71.jar \
        run /home/mbrown/GitHub/VirusFinder2_VERSE/WDL/VirusFinder2.wdl \
        -i simulation_inputs_virusfinder.json
