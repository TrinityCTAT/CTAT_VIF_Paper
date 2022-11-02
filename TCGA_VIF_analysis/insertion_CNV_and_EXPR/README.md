

Simplify insertions, make non-redundant:

- first restrict to highest evidence insertion per virus breakpoint group.
- then each sample and within 1Mb region, select the single highest evidence insertion.

```
    ./simplify_best_insertion_per_sample_region.py  --insertions_tsv ../filtered_insertions.tsv --output filtered_insertions.nr.tsv
```


