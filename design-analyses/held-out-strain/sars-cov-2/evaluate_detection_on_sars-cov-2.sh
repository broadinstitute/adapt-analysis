#!/bin/bash

# Evaluate the design against SARS-CoV-2 genomes up to 2020-11-12.

# Load ADAPT conda environment and set some environment variables
source ~/anaconda3/etc/profile.d/conda.sh
conda activate adapt
PREDICTIVE_MODELS_PATH="/home/hayden/adapt/models"
ARG_PM="3"

GENOMES_PATH="./input/gisaid_msa_2020-11-12/msa_1112.fasta.gz"

function evaluate {
    # Args:
    #   1: name of design

    # Determine active fraction
    analyze_coverage.py designs/${1}.tsv.0 $GENOMES_PATH -pm $ARG_PM --write-frac-bound evaluations/${1}.sars-cov-2.frac-bound.active.tsv --verbose --predict-activity-model-path $PREDICTIVE_MODELS_PATH/classify/model-51373185 $PREDICTIVE_MODELS_PATH/regress/model-f8b6fd5d

    # Determine highly active fraction
    analyze_coverage.py designs/${1}.tsv.0 $GENOMES_PATH -pm $ARG_PM --write-frac-bound evaluations/${1}.sars-cov-2.frac-bound.highly-active.tsv --verbose --predict-activity-model-path $PREDICTIVE_MODELS_PATH/classify/model-51373185 $PREDICTIVE_MODELS_PATH/regress/model-f8b6fd5d --predict-activity-require-highly-active

    # Determine mean activity of each guide set
    analyze_coverage.py designs/${1}.tsv.0 $GENOMES_PATH -pm $ARG_PM --write-mean-activity-of-guides evaluations/${1}.sars-cov-2.mean-guide-activities.tsv --verbose --predict-activity-model-path $PREDICTIVE_MODELS_PATH/classify/model-51373185 $PREDICTIVE_MODELS_PATH/regress/model-f8b6fd5d
}

evaluate "design"
evaluate "design.downsampled-sars-cov-1"
