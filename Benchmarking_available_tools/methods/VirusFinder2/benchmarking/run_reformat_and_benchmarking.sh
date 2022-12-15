#!/bin/bash

set -ex


./scripts/run_benchmarking_VF2.py \
    --Directory `pwd`/data \
    --insertion_truth_set `pwd`/insertion_truth_set.tsv


