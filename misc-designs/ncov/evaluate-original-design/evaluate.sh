#!/bin/bash

# Evaluate the original, Jan 22 2020 assay design against
# the latest genomes.
# ADAPT must be loaded (conda env) before running.

ADAPT_PATH="$HOME/adapt"
GENOMES_PATH="./data/gisaid_msa_2021-03-04/msa_0304.fasta.gz"

# Determine active fraction
analyze_coverage.py sars-cov-2.v1-experimentally-tested.tsv $GENOMES_PATH -pm 2 --write-frac-bound frac-bound.active.tsv --verbose --predict-activity-model-path $ADAPT_PATH/models/classify/model-51373185 $ADAPT_PATH/models/regress/model-f8b6fd5d

# Determine highly active fraction
analyze_coverage.py sars-cov-2.v1-experimentally-tested.tsv $GENOMES_PATH -pm 2 --write-frac-bound frac-bound.highly-active.tsv --verbose --predict-activity-model-path $ADAPT_PATH/models/classify/model-51373185 $ADAPT_PATH/models/regress/model-f8b6fd5d --predict-activity-require-highly-active
