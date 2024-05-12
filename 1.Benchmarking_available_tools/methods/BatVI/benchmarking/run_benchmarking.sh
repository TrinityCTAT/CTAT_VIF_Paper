#!/bin/bash

set -ex

./scripts/run_benchmarking_batvi.py \
	--Directory `pwd`/data \
	--insertion_truth_set `pwd`/insertion_truth_set.tsv

