#!/bin/bash

# Compile evaluation files into one.

OUT_FN="evaluations/evaluations-compiled.tsv"

# Write header
echo -e "design\tevaluated_against\tdesign_id\tstat\tvalue" > $OUT_FN

function copy_cols {
    # Args:
    #   1: prefix of filename
    #   2: design name
    #   3: evaluated against

    in_prefix="evaluations/$1"

    # Frac bound (active)
    cat "${in_prefix}.frac-bound.active.tsv" | tail -n +2 | awk -F'\t' -v design="$2" -v against="$3" '{print design"\t"against"\t"$1"\tfrac-bound-active\t"$5}' >> $OUT_FN

    # Frac bound (highly active)
    cat "${in_prefix}.frac-bound.highly-active.tsv" | tail -n +2 | awk -F'\t' -v design="$2" -v against="$3" '{print design"\t"against"\t"$1"\tfrac-bound-highly-active\t"$5}' >> $OUT_FN

    # Mean guide activities
    cat "${in_prefix}.mean-guide-activities.tsv" | tail -n +2 | awk -F'\t' -v design="$2" -v against="$3" '{print design"\t"against"\t"$1"\tmean-guide-activity\t"$3}' >> $OUT_FN
}

# Designs with all data
copy_cols "design.sars-cov-2" "all" "sars-cov-2"
copy_cols "design.design-input" "all" "input"

# Designs with downsampled SARS-CoV-1 data
copy_cols "design.downsampled-sars-cov-1.sars-cov-2" "downsampled-sars-cov-1" "sars-cov-2"
copy_cols "design.downsampled-sars-cov-1.design-input" "downsampled-sars-cov-1" "input"
