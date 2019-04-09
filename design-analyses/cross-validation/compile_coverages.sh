#!/bin/bash

# Compile a distribution of coverages for each taxonomy.
#
# Author: Hayden Metsky <hayden@mit.edu>


# Let the coverage of a design be the mean of the top TOPN targets
# in the design
TOPN=10

for dir in $(ls -1 . | grep 'tax-'); do
    design_name=$(echo "$dir" | sed -s 's/tax-//')
    echo "Compiling coverages for $design_name"

    # Start a file with the distribution of coverage on the test set
    echo -n "" > tax-${design_name}/coverage-against-test.distribution.txt

    # Compute a coverage value for each resampling
    for design_file in $(ls -1 tax-${design_name}/designs/coverages/design-*.txt); do
        topn_mean=$(cat $design_file | tail -n +2 | awk -F '\t' -v topn="$TOPN" '$1<=topn {s+=$5} END {print s/topn}')
        echo "$topn_mean" >> tax-${design_name}/coverage-against-test.distribution.txt
    done
done
