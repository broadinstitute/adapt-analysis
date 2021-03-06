#!/bin/bash

# Evaluate dispersion from designs.
#
# Author: Hayden Metsky <hayden@mit.edu>


for dir in $(ls -1 . | grep 'tax-'); do
    design_name=$(echo "$dir" | sed -s 's/tax-//')
    echo "Evaluating dispersion for $design_name"

    # Summarize dispersion when input is resampled sequences
    python ../scripts/compute_stats_from_design.py dispersion tax-${design_name}/designs/resampled/design-*.tsv.0 --num-targets 5 --pairwise-jaccard-distribution-out tax-${design_name}/dispersion.resampled.distribution.txt --pairwise-jaccard-distribution-with-loose-equality-out tax-${design_name}/dispersion.resampled.loose-equality.distribution.txt > tax-${design_name}/dispersion.resampled.summary.txt

    # Summarize dispersion when input is all (non-resampled) sequences
    python ../scripts/compute_stats_from_design.py dispersion tax-${design_name}/designs/non-resampled/design-*.tsv.0 --num-targets 5 --pairwise-jaccard-distribution-out tax-${design_name}/dispersion.non-resampled.distribution.txt --pairwise-jaccard-distribution-with-loose-equality-out tax-${design_name}/dispersion.non-resampled.loose-equality.distribution.txt > tax-${design_name}/dispersion.non-resampled.summary.txt
done
