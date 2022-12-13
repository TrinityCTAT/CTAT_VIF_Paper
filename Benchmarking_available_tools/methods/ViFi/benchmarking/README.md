# format conversion
    
./scripts/reformatViFi.py --inputs ./data/Test.clusters.txt --sampleID subset.insertion_seqs.fa  --output `pwd`/benchmarking_out

# benchmark

./scripts/Benchmarking.py  --truth_insertions insertion_truth_set.tsv  --predicted_insertions `pwd`/benchmarking_out/ViFi_Output_reformated.txt --output `pwd`/benchmarking_out
