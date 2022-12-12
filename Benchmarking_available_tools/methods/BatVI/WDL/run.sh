#!/bin/bash

set -e


sudo rm -r cromwell*

#java -jar /home/mbrown/Tools/cromwell-71.jar \
#	run /home/mbrown/GitHub/BatVI/WDL/BatVI.wdl \
#	-i inputs.json


java -jar /home/mbrown/Tools/cromwell-71.jar \
        run BatVI.wdl \
        -i inputs.json
