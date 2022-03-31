#!/bin/bash

# adapt conda environment should be activated before running
# Args:
#   1: name of taxon (e.g., 'lasv-L')

python $HOME/adapt/bin/analyze_coverage.py data/$1/assay-options.tsv data/$1/aln.fasta --write-per-seq out/$1/per-seq-predictions.tsv -pm 5 --predict-cas13a-activity-model latest latest
