#!/bin/bash

# Design specific panel.

# Set input/output
IN_TSV="input/coronaviruses.tsv"
OUT_TSV_DIR="out/designs"

MANUAL_SEQ_INPUT="input/manual-fasta-specify.tsv"
ONLY_DESIGN_NCOV="input/only-design.ncov.tsv"
ONLY_DESIGN_OTHER_HUMAN="input/only-design.other-human.tsv"

# Set variables for design
PREP_MEMOIZE_DIR="/ebs/dgd-analysis/prep-memoize-dir"
MAFFT_PATH="/home/hayden/viral-ngs/viral-ngs-etc/conda-env/bin/mafft"
CLUSTER_THRESHOLD=0.20
ARG_GL="28"
ARG_GM="1"
ARG_GP="0.99"
ARG_PL="30"
ARG_PM="1"
ARG_PP="0.99"
ARG_IDM="5"
ARG_IDFRAC="0"
ARG_MAXPRIMERSATSITE="5"
ARG_MAXTARGETLENGTH="250"
ARG_COSTFNWEIGHTS="0.6667 0.2222 0.1111"
ARG_BESTNTARGETS="200"

# Set more strict variables for nCoV
ARG_GM_NCOV="0"
ARG_PM_NCOV="0"

# Set API key
# Set as environment variable
NCBI_API_KEY=""

# Set a predictive model
# This is from commit da10963 of my adapt-seq-design repo
PREDICTIVE_MODEL="/home/hayden/adapt-seq-design/models/predictor_exp-and-pos_regress-on-active/model-8f534a8c"

# Make the memoize directory
mkdir -p $PREP_MEMOIZE_DIR

# Set tmp directory
export TMPDIR="/tmp"

# Activate adapt conda environment
source ~/anaconda3/etc/profile.d/conda.sh
conda activate /ebs/dgd-analysis/tools/envs/adapt-prod

mkdir -p out/designs

# Split the design into two parts, allowing more strict parameters for nCoV design

# Run the design, only for ncov
/usr/bin/time -f "mem=%K RSS=%M elapsed=%e cpu.sys=%S .user=%U" design.py complete-targets auto-from-file $IN_TSV $OUT_TSV_DIR -gl $ARG_GL -gm $ARG_GM_NCOV -gp $ARG_GP -pl $ARG_PL -pm $ARG_PM_NCOV -pp $ARG_PP --id-m $ARG_IDM --id-frac $ARG_IDFRAC --id-method shard --require-flanking3 H --max-primers-at-site $ARG_MAXPRIMERSATSITE --primer-gc-content-bounds 0.4 0.6 --max-target-length $ARG_MAXTARGETLENGTH --cost-fn-weights $ARG_COSTFNWEIGHTS --best-n-targets $ARG_BESTNTARGETS --predict-activity-model-path $PREDICTIVE_MODEL --mafft-path $MAFFT_PATH --prep-memoize-dir $PREP_MEMOIZE_DIR --cluster-threshold $CLUSTER_THRESHOLD --use-fasta $MANUAL_SEQ_INPUT --only-design-for $ONLY_DESIGN_NCOV --write-input-seqs --write-input-aln --ncbi-api-key $NCBI_API_KEY --verbose &> out/design.ncov.out

# Run the design, for other human coronaviruses
/usr/bin/time -f "mem=%K RSS=%M elapsed=%e cpu.sys=%S .user=%U" design.py complete-targets auto-from-file $IN_TSV $OUT_TSV_DIR -gl $ARG_GL -gm $ARG_GM -gp $ARG_GP -pl $ARG_PL -pm $ARG_PM -pp $ARG_PP --id-m $ARG_IDM --id-frac $ARG_IDFRAC --id-method shard --require-flanking3 H --max-primers-at-site $ARG_MAXPRIMERSATSITE --primer-gc-content-bounds 0.4 0.6 --max-target-length $ARG_MAXTARGETLENGTH --cost-fn-weights $ARG_COSTFNWEIGHTS --best-n-targets $ARG_BESTNTARGETS --predict-activity-model-path $PREDICTIVE_MODEL --mafft-path $MAFFT_PATH --prep-memoize-dir $PREP_MEMOIZE_DIR --cluster-threshold $CLUSTER_THRESHOLD --use-fasta $MANUAL_SEQ_INPUT --only-design-for $ONLY_DESIGN_OTHER_HUMAN --write-input-seqs --write-input-aln --ncbi-api-key $NCBI_API_KEY --verbose &> out/design.other-human.out

# gzip the stdout/stderr
gzip out/design.ncov.out
gzip out/design.other-human.out
