#!/bin/bash

# Set variables for design
PREP_MEMOIZE_DIR="/ebs/dgd-analysis/prep-memoize-dir"
MAFFT_PATH="/home/hayden/viral-ngs/viral-ngs-etc/conda-env/bin/mafft"
CLUSTER_THRESHOLD=1.0
ARG_GL="28"
ARG_GM="1"
ARG_GP="0.99"
ARG_PL="30"
ARG_PM="1"
ARG_PP="0.99"
ARG_MAXPRIMERSATSITE="5"
ARG_MAXTARGETLENGTH="250"
ARG_COSTFNWEIGHTS="0.6667 0.2222 0.1111"
ARG_BESTNTARGETS="20"

# Set API key
NCBI_API_KEY="321e5004f2502f604ee2e8858a22b993d608"

# Set a predictive model
# This is from commit da10963 of my adapt-seq-design repo
PREDICTIVE_MODEL="/ebs/dgd-analysis/tools/adapt-seq-design-prod/models/predictor_exp-and-pos_regress-on-active/model-8f534a8c"

# Make the memoize directory
mkdir -p $PREP_MEMOIZE_DIR

# Set tmp directory
export TMPDIR="/tmp"

# Activate adapt conda environment
source ~/anaconda3/etc/profile.d/conda.sh
conda activate /ebs/dgd-analysis/tools/envs/adapt-prod

OUT_TSV_DIR="out/designs"
mkdir -p $OUT_TSV_DIR
