This designs for the 6 human-infecting coronavirus species with >= 10 genomes (including nCoV), such that each design is specific against the other 47 coronavirus species.
This uses `--id-m 4 --id-frac 0` to enforce specificity.

`out-lsh` contains output using the LSH near-neighbor method to query for specificity. `out` contains the method sharding k-mers across tries.
