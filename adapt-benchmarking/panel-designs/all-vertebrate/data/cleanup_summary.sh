#!/bin/sh

IN="summary.one-segment-per-taxid.tsv"
OUT="design-results.tsv"

# Write header
echo "family\tgenus\tspecies\ttaxid\tsegment\texperiment\tresult_type\tcluster\ttimestamp\tlast_changed_timestamp\tnum_input_seqs\tnum_curated_seqs\tnum_clusters\trss\telapsed_time\tmean_cluster_objective_value\tmean_cluster_target_len\tmean_cluster_num_primers5\tmean_cluster_num_primers3\tmean_cluster_num_guides\tmean_cluster_guide_set_frac_bound\tmean_cluster_guide_set_expected_activity\tmean_cluster_guide_set_median_activity\tmean_cluster_guide_set_5th_pctile_activity" > $OUT

# Copy rows on taxons (skip the per-cluster rows)
cat $IN | awk -F'\t' '$7=="taxon" {print $0}' >> $OUT
