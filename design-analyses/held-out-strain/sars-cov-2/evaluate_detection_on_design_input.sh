#!/bin/bash

# Evaluate the design against the input used to make the design
# (i.e., how well would the designs do in the usual case where
# they are used against what they are designed against).

# Load ADAPT conda environment and set some environment variables
source ~/anaconda3/etc/profile.d/conda.sh
conda activate adapt
PREDICTIVE_MODELS_PATH="/home/hayden/adapt/models"
ARG_PM="3"

GENOMES_PATH_ALL="./input/sars-related-cov_upto_2019-11-30.fasta.gz"
GENOMES_PATH_DOWNSAMPLED_SARS_COV_1="./input/sars-related-cov_upto_2019-11-30.single-sars-cov-1.fasta.gz"

function evaluate {
    # Args:
    #   1: name of design
    #   2: path to sequences to evaluate against

    # Determine active fraction
    analyze_coverage.py designs/${1}.tsv.0 $2 -pm $ARG_PM --write-frac-bound evaluations/${1}.design-input.frac-bound.active.tsv --verbose --predict-activity-model-path $PREDICTIVE_MODELS_PATH/classify/model-51373185 $PREDICTIVE_MODELS_PATH/regress/model-f8b6fd5d

    # Determine highly active fraction
    analyze_coverage.py designs/${1}.tsv.0 $2 -pm $ARG_PM --write-frac-bound evaluations/${1}.design-input.frac-bound.highly-active.tsv --verbose --predict-activity-model-path $PREDICTIVE_MODELS_PATH/classify/model-51373185 $PREDICTIVE_MODELS_PATH/regress/model-f8b6fd5d --predict-activity-require-highly-active

    # Determine mean activity of each guide set
    analyze_coverage.py designs/${1}.tsv.0 $2 -pm $ARG_PM --write-mean-activity-of-guides evaluations/${1}.design-input.mean-guide-activities.tsv --verbose --predict-activity-model-path $PREDICTIVE_MODELS_PATH/classify/model-51373185 $PREDICTIVE_MODELS_PATH/regress/model-f8b6fd5d
}

evaluate "design" $GENOMES_PATH_ALL
evaluate "design.downsampled-sars-cov-1" $GENOMES_PATH_DOWNSAMPLED_SARS_COV_1
