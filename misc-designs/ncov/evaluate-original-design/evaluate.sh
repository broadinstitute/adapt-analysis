#!/bin/bash

# Evaluate the original, Jan 22 2020 assay design against
# the latest genomes.
# ADAPT must be loaded (conda env) before running.

ADAPT_PATH="$HOME/adapt"
GENOMES_PATH="./data/gisaid_msa_2021-12-12/msa_1212.fasta.gz"

# Determine active fraction
analyze_coverage.py sars-cov-2.v1-experimentally-tested.tsv $GENOMES_PATH -pm 2 --write-frac-bound frac-bound.active.tsv --verbose --predict-activity-model-path $ADAPT_PATH/adapt/models/classify/cas13a/v1_0 $ADAPT_PATH/adapt/models/regress/cas13a/v1_0

# Determine highly active fraction
analyze_coverage.py sars-cov-2.v1-experimentally-tested.tsv $GENOMES_PATH -pm 2 --write-frac-bound frac-bound.highly-active.tsv --verbose --predict-activity-model-path $ADAPT_PATH/adapt/models/classify/cas13a/v1_0 $ADAPT_PATH/adapt/models/regress/cas13a/v1_0 --predict-activity-require-highly-active
