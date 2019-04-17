#!/bin/bash

# Compile a distribution of coverages for each taxonomy.
#
# Author: Hayden Metsky <hayden@mit.edu>


# Let the coverage of a design be the mean of the top TOPN targets
# in the design
TOPN=20

for dir in $(ls -1 . | grep 'tax-'); do
    design_name=$(echo "$dir" | sed -s 's/tax-//')
    echo "Compiling coverages for $design_name"

    # Start a file with the distribution of coverage on the test set
    echo -n "" > tax-${design_name}/coverage-against-test.distribution.txt

    # Compute a coverage value for each resampling
    for design_file in $(ls -1 tax-${design_name}/designs/coverages/design-*.txt); do
        # Last awk command takes median; everything before extracts/sorts the top TOPN coverage values
        #topn_median=$(cat $design_file | tail -n +2 | awk -F '\t' -v topn="$TOPN" '$1<=topn {print $5}' | sort -g | awk '{ covg[NR] = $1 } END { if (NR % 2) { print covg[(NR + 1) / 2] } else { print (covg[(NR / 2)] + covg[(NR / 2) + 1]) / 2.0 } }')
        topn_mean=$(cat $design_file | tail -n +2 | awk -F '\t' -v topn="$TOPN" '$1<=topn {s+=$5; n+=1} END {print s/n}')
        echo "$topn_mean" >> tax-${design_name}/coverage-against-test.distribution.txt
    done
done
