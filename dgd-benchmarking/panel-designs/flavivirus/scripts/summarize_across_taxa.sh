#!/bin/bash

# Summarize designs across taxonomies.

# This is a far simpler version of ~/misc-repos/continuous-dgd/summarize_across_taxa.sh,
# since no species are segmented and we only have 1 cluster per species.

TAXONOMIES_FILE="flaviviruses.clade.tsv"

# Write header
echo -e "name\ttaxid\tsegment\tcost\ttarget_len\tnum_primers5\tnum_primers3\tnum_guides\tnonspecific_cost\tnonspecific_target_len\tnonspecific_num_primers5\tnonspecific_num_primers3\tnonspecific_num_guides" > out/summary.tsv

while read -r taxonomy; do
    name=$(echo "$taxonomy" | awk -F '\t' '{print $1}')
    taxid=$(echo "$taxonomy" | awk -F'\t' '{print $2}')
    segment=$(echo "$taxonomy" | awk -F'\t' '{print $3}')
    refaccs=$(echo "$taxonomy" | awk -F'\t' '{print $4}')

    # Find stats for this design, for the best target (first one)
    cost=$(cat out/designs/${name}.0.tsv | tail -n +2 | head -n 1 | awk -F'\t' '{print $1}')
    target_len=$(cat out/designs/${name}.0.tsv | tail -n +2 | head -n 1 | awk -F'\t' '{print $4}')
    num_primers5=$(cat out/designs/${name}.0.tsv | tail -n +2 | head -n 1 | awk -F'\t' '{print $6}')
    num_primers3=$(cat out/designs/${name}.0.tsv | tail -n +2 | head -n 1 | awk -F'\t' '{print $10}')
    num_guides=$(cat out/designs/${name}.0.tsv | tail -n +2 | head -n 1 | awk -F'\t' '{print $13}')

    # Write stats
    echo -ne "$name\t$taxid\t$segment\t$cost\t$target_len\t$num_primers5\t$num_primers3\t$num_guides\t" >> out/summary.tsv

    # Get the summary line from the larger design, without
    # specificity (if the taxid was in that design -- it may have
    # not if it has < 10 seqs)
    nonspecific=$(cat ../all-viral-with-ge10-seqs/data/design-results.tsv | awk -F'\t' -v taxid="$taxid" '$4==taxid {print $0}')
    if [ ! -z "$nonspecific" ]; then
        # Only do this for species in that large design
        nonspecific_cost=$(echo "$nonspecific" | awk -F'\t' '{print $16}')
        nonspecific_target_len=$(echo "$nonspecific" | awk -F'\t' '{print $17}')
        nonspecific_num_primers5=$(echo "$nonspecific" | awk -F'\t' '{print $18}')
        nonspecific_num_primers3=$(echo "$nonspecific" | awk -F'\t' '{print $19}')
        nonspecific_num_guides=$(echo "$nonspecific" | awk -F'\t' '{print $20}')

        # Write stats
        echo -e "$nonspecific_cost\t$nonspecific_target_len\t$nonspecific_num_primers5\t$nonspecific_num_primers3\t$nonspecific_num_guides" >> out/summary.tsv
    else
        # Write 'NA' for these fields
        echo -e "NA\tNA\tNA\tNA\tNA" >> out/summary.tsv
    fi
done < <(tail -n +2 $TAXONOMIES_FILE)
