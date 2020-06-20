#!/bin/bash

# Compute a few simple summary statistics on the designs.

IN="data/design-results.tsv"
HIGHLY_ACTIVE_THRES="2.7198637"

div() { awk -v a="$1" -v b="$2" "BEGIN { print a/b }"; }

# Only use the 'specific_max-activity' designs
data=$(mktemp)
cat $IN | awk -F'\t' '$6=="specific_max-activity" {print $0}' > $data

# Compute stats on number of designs/species
# The number of designs and number of species should be the same, but
# output both to check
# (we should have already filtered to 1 segment per species)
num_designs=$(cat $data | wc -l)
num_species=$(cat $data | cut -f4 | sort | uniq | wc -l)
echo "Number of designs: $num_designs"
echo "Number of species: $num_species"
echo "  (These numbers should be the same)"
echo ""

# Compute stats on input size
num_input_seqs=$(cat $data | cut -f11 | awk '{s += $1} END {print s}')
echo "Total number of input sequences: $num_input_seqs"
echo "  (Since 1 segment per species, this is roughly the total number of input genomes)"
echo ""

# Compute stats on predicted activity
num_designs_median_high=$(cat $data | cut -f23 | awk -v t="$HIGHLY_ACTIVE_THRES" '$1 >= t {s += 1} END {print s}')
frac_designs_median_high=$(div "$num_designs_median_high" "$num_designs")
num_designs_5th_pctile_high=$(cat $data | cut -f24 | awk -v t="$HIGHLY_ACTIVE_THRES" '$1 >= t {s += 1} END {print s}')
frac_designs_5th_pctile_high=$(div "$num_designs_5th_pctile_high" "$num_designs")
echo "Fraction of designs where median activity is high: $frac_designs_median_high"
echo "Fraction of designs where 5th pctile activity is high: $frac_designs_5th_pctile_high"
echo ""

# Compute stats on runtime
num_designs_under_24hrs=$(cat $data | cut -f15 | awk '$1 <= 60*60*24 {s += 1} END {print s}')
num_designs_under_4hrs=$(cat $data | cut -f15 | awk '$1 <= 60*60*4 {s += 1} END {print s}')
echo "Number of designs completed in <= 24 hrs: $num_designs_under_24hrs"
echo "Number of designs completed in <= 4 hrs: $num_designs_under_4hrs"
echo ""

# Compute stats on curation
num_designs_removed_50pct_seqs=$(cat $data | cut -f11,12 | awk '$2/$1 < 0.5 {s += 1} END {print s}')
echo "Number of designs where >50% of sequences are removed during curation: $num_designs_removed_50pct_seqs"

rm $data
