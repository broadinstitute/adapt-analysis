#!/bin/sh

IN="summary.one-segment-per-taxid.tsv"
OUT="design-results.tsv"

# Write header
echo "family\tgenus\tspecies\ttaxid\tsegment\tnum_seqs\tresult_type\tcluster\ttimestamp\tlast_changed_timestamp\tnum_input_seqs\tnum_curated_seqs\tnum_clusters\trss\telapsed_time\tmean_cluster_cost\tmean_cluster_target_len\tmean_cluster_num_primers5\tmean_cluster_num_primers3\tmean_cluster_num_guides" > $OUT

# Copy rows on taxons (skip the per-cluster rows)
cat $IN | awk -F'\t' '$7=="taxon" {print $0}' >> $OUT
