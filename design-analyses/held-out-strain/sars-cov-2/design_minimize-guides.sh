#!/bin/bash

# Design for SARS-related CoV species using all sequences
# up to 2019-11-30 (i.e., before SARS-CoV-2).
# Instead of maximizing activity, this minimizes the assay
# complexity subject to constraints on coverage.

# Load run script from adapt-designs repo
source ~/misc-repos/adapt-designs/scripts/run-adapt/run_common.sh

# Make a file to tell ADAPT to use a fasta containing
# the sequences, rather than downloading the most recent
usefasta=$(mktemp)
echo -e "694009\tNone\tinput/sars-related-cov_upto_2019-11-30.fasta.gz" > $usefasta

taxid="694009"
segment="None"
ref_acc="NC_004718"

# Change some parameters from their default values
CLUSTER_THRESHOLD="0.30"
ARG_PP="0.99"
ARG_GP="0.99"
ARG_MAXPRIMERSATSITE="10"

OUT_DIR="designs"
mkdir -p $OUT_DIR

specificity_taxa="input/specificity_taxa.tsv"
ARG_SPECIFICITY="--id-m $ARG_IDM --id-frac $ARG_IDFRAC --id-method shard --specific-against-taxa $specificity_taxa"

ARG_OBJ="--obj minimize-guides -gp $ARG_GP"

# Run the design
run-adapt design_minimize-guides complete-targets auto-from-args $taxid "$segment" $ref_acc $OUT_DIR/design.minimize-guides.tsv -gl $ARG_GL -pl $ARG_PL -pm $ARG_PM -pp $ARG_PP --primer-gc-content-bounds $ARG_PRIMER_GC_LO $ARG_PRIMER_GC_HI --max-primers-at-site $ARG_MAXPRIMERSATSITE --max-target-length $ARG_MAXTARGETLENGTH $ARG_OBJ --obj-fn-weights $ARG_OBJFNWEIGHTS --best-n-targets $ARG_BESTNTARGETS $ARG_SPECIFICITY --predict-activity-model-path $ARG_PREDICTIVE_MODELS --mafft-path $MAFFT_PATH --prep-memoize-dir $PREP_MEMOIZE_DIR --ncbi-api-key $NCBI_API_KEY --cluster-threshold $CLUSTER_THRESHOLD --use-fasta $usefasta --seed 1 --verbose

# Remove the tmp file
rm $usefasta