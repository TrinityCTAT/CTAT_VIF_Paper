#!/bin/bash

set -ex


docker build --build-arg CACHEBUST=$(date +%s) -f Dockerfile -t trinityctat/nf_vif:devel .

