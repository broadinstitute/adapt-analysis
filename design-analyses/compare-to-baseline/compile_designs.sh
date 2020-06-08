#!/bin/bash

# Compile a distribution of designs (real and naive) for each taxonomy.
#
# Each design is from a sliding window. The real min-guides designs achieve a desired coverage
# (e.g., 99%) of the targets with one or more guides. The real max-activity designs maximize
# coverage (expected activity, where activity is 1 or 0) for different constraints on the
# number of guides. The naive designs use only one guide per window, and report on the coverage
# obtained with the guide.
# The designs are produced multiple times, each with a separate resampling of the
# input; this way, we can produce a distribution (e.g., of coverage obtained by
# the naive design) for each window.
#
# Author: Hayden Metsky <hayden@mit.edu>


function compile_for_taxid() {
    taxid="$1"
    segment="$2"
    tax_dir="tax-${taxid}_${segment}"
    echo "Compiling designs for $taxid (segment: $segment)"

    # Start a file with the distribution of designs for the real designs that maximize activity (coverage)
    echo -n "" > $tax_dir/real-designs.max-activity.tsv
    # Add the header to the file, and add column at end giving guide constraint (hgc)
    head -n 1 $tax_dir/designs/design-1.real-design.maximize-activity.hgc-1.tsv | awk '{print $0"\thgc"}' >> $tax_dir/real-designs.max-activity.tsv
    # Add to the file the designs (for each window) from each resampling
    for hgc in $(seq 1 5); do
        for design_file in $(ls -1 $tax_dir/designs/design-*.real-design.maximize-activity.hgc-${hgc}.tsv); do
            # Leave out the header when reading the file, and add column with hgc
            tail -n +2 ${design_file} | awk -v hgc="$hgc" '{print $0"\t"hgc}' >> $tax_dir/real-designs.max-activity.tsv
        done
    done

    # Start a file with the distribution of designs for the real designs that minimize guides
    echo -n "" > $tax_dir/real-designs.min-guides.tsv
    # Add the header to the file
    head -n 1 $tax_dir/designs/design-1.real-design.minimize-guides.tsv >> $tax_dir/real-designs.min-guides.tsv
    # Add to the file the designs (for each window) from each resampling
    for design_file in $(ls -1 $tax_dir/designs/design-*.real-design.minimize-guides.tsv); do
        # Leave out the header when reading the file
        tail -n +2 ${design_file} >> $tax_dir/real-designs.min-guides.tsv
    done

    # Start a file with the distribution of designs for the naive designs
    echo -n "" > $tax_dir/naive-designs.tsv
    # Add the header to the file
    head -n 1 $tax_dir/designs/design-1.naive-design.tsv >> $tax_dir/naive-designs.tsv
    # Add to the file the designs (for each window) from each resampling
    for design_file in $(ls -1 $tax_dir/designs/design-*.naive-design.tsv); do
        # Leave out the header when reading the file
        tail -n +2 ${design_file} >> $tax_dir/naive-designs.tsv
    done

    # gzip the files
    gzip -f $tax_dir/real-designs.max-activity.tsv
    gzip -f $tax_dir/real-designs.min-guides.tsv
    gzip -f $tax_dir/naive-designs.tsv
}

compile_for_taxid "64320" "None"
compile_for_taxid "11620" "S"
compile_for_taxid "11620" "L"
compile_for_taxid "121791" "None"
compile_for_taxid "11676" "None"
compile_for_taxid "11103" "None"
compile_for_taxid "11320" "2"
compile_for_taxid "147711" "None"
compile_for_taxid "147712" "None"
compile_for_taxid "463676" "None"
compile_for_taxid "138948" "None"
compile_for_taxid "138949" "None"
compile_for_taxid "138950" "None"
compile_for_taxid "138951" "None"
