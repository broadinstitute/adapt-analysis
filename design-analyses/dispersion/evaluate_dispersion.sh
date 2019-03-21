#!/bin/bash

# Evaluate dispersion from designs.
#
# Author: Hayden Metsky <hayden@mit.edu>


for dir in $(ls -1 . | grep 'tax-'); do
    design_name=$(echo "$dir" | sed -s 's/tax-//')
    echo "Evaluating dispersion for $design_name"

    # Summarize dispersion
    python ../scripts/compute_stats_from_design.py dispersion tax-${design_name}/designs/design-*.tsv.0 > tax-${design_name}/dispersion.summary.txt
done
