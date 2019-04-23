#!/bin/bash

# Compile a distribution of designs (real and naive) for each taxonomy.
#
# Each design is from a sliding window. The real designs achieve a desired coverage
# (e.g., 99%) of the targets with one or more guides. Th naive designs use only
# one guide per window, and report on the coverage obtained with the guide.
# The designs are produced multiple times, each with a separate resampling of the
# input; this way, we can produce a distribution (e.g., of coverage obtained by
# the naive design) for each window.
#
# Author: Hayden Metsky <hayden@mit.edu>


for dir in $(ls -1 . | grep 'tax-'); do
    tax_name=$(echo "$dir" | sed -s 's/tax-//')
    echo "Compiling designs for $tax_name"

    # Start a file with the distribution of designs for the real designs
    echo -n "" > tax-${tax_name}/real-designs.tsv
    # Add the header to the file
    head -n 1 tax-${tax_name}/designs/design-1.real-design.tsv.0 >> tax-${tax_name}/real-designs.tsv
    # Add to the file the designs (for each window) from each resampling
    for design_file in $(ls -1 tax-${tax_name}/designs/design-*.real-design.tsv.0); do
        # Leave out the header when reading the file
        tail -n +2 ${design_file} >> tax-${tax_name}/real-designs.tsv
    done

    # Start a file with the distribution of designs for the naive designs
    echo -n "" > tax-${tax_name}/naive-designs.tsv
    # Add the header to the file
    head -n 1 tax-${tax_name}/designs/design-1.naive-design.tsv >> tax-${tax_name}/naive-designs.tsv
    # Add to the file the designs (for each window) from each resampling
    for design_file in $(ls -1 tax-${tax_name}/designs/design-*.naive-design.tsv); do
        # Leave out the header when reading the file
        tail -n +2 ${design_file} >> tax-${tax_name}/naive-designs.tsv
    done
done
