#!/bin/bash

# Summarize designs across taxonomies.

# This is a far simpler version of ~/misc-repos/continuous-dgd/summarize_across_taxa.sh,
# since no species are segmented and we only have 1 cluster per species.

TAXONOMIES_FILE="flaviviruses.clade.tsv"

# Write header
echo -e "name\ttaxid\tsegment\tcost\ttarget_len\tnum_primers5\tnum_primers3\tnum_guides\tnonspecific_cost\tnonspecific_target_len\tnonspecific_num_primers5\tnonspecific_num_primers3\tnonspecific_num_guides" > summary/summary.tsv

while read -r taxonomy; do
    name=$(echo "$taxonomy" | awk -F '\t' '{print $1}')
    taxid=$(echo "$taxonomy" | awk -F'\t' '{print $2}')
    segment=$(echo "$taxonomy" | awk -F'\t' '{print $3}')
    refaccs=$(echo "$taxonomy" | awk -F'\t' '{print $4}')

    # For all species, the parameters in out/ match what was used
    #  for the non-specific design, so pull out these so we can
    #  compare with that (except for some parameters that don't matter --
    #  such as the number N of best designs to find, or the max number of
    #  primers at a site (as long as all designs falls under this))
    # EXCEPT for dengue: for this, the parameters in the non-specific design
    #  were a little different (namely, ARG_PM=2 rather than ARG_PM=3), so
    #  for dengue pull out the specific design using these parameters (it
    #  shouldn't matter that there's a difference in the value N for these
    #  two runs, since we always just parse the best design -- so all that
    #  was really needed is N=1)
    if [ $taxid == "12637" ]; then
        out_dir="out_PM2_N5"
    else
        out_dir="out"
    fi

    # Find stats for this design, for the best target (first one)
    cost=$(cat ${out_dir}/designs/${name}.0.tsv | tail -n +2 | head -n 1 | awk -F'\t' '{print $1}')
    target_len=$(cat ${out_dir}/designs/${name}.0.tsv | tail -n +2 | head -n 1 | awk -F'\t' '{print $4}')
    num_primers5=$(cat ${out_dir}/designs/${name}.0.tsv | tail -n +2 | head -n 1 | awk -F'\t' '{print $6}')
    num_primers3=$(cat ${out_dir}/designs/${name}.0.tsv | tail -n +2 | head -n 1 | awk -F'\t' '{print $10}')
    num_guides=$(cat ${out_dir}/designs/${name}.0.tsv | tail -n +2 | head -n 1 | awk -F'\t' '{print $13}')

    # Write stats
    echo -ne "$name\t$taxid\t$segment\t$cost\t$target_len\t$num_primers5\t$num_primers3\t$num_guides\t" >> summary/summary.tsv

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
        echo -e "$nonspecific_cost\t$nonspecific_target_len\t$nonspecific_num_primers5\t$nonspecific_num_primers3\t$nonspecific_num_guides" >> summary/summary.tsv
    else
        # Write 'NA' for these fields
        echo -e "NA\tNA\tNA\tNA\tNA" >> summary/summary.tsv
    fi
done < <(tail -n +2 $TAXONOMIES_FILE)
