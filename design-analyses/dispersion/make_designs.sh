#!/bin/bash

# Make designs to use for evaluating dispersion.
#
# Author: Hayden Metsky <hayden@mit.edu>


# Load environment, and variables, for ADAPT
source ~/misc-repos/adapt-designs/scripts/run-adapt/custom-env/load_custom_env.sh

# Set variables for measuring dispersion
NUM_DESIGNS=20
NJOBS=10

# Set variables for design
CLUSTER_THRESHOLD=1.0   # Use high value to obtain a single cluster
ARG_GL="28"
ARG_PL="30"
ARG_PM="3"
ARG_PP="0.98"
ARG_PRIMER_GC_LO="0.35"
ARG_PRIMER_GC_HI="0.65"
ARG_SOFTGUIDECONSTRAINT="1"
ARG_HARDGUIDECONSTRAINT="5"
ARG_PENALTYSTRENGTH="0.25"
ARG_MAXIMIZATIONALGORITHM="random-greedy"
ARG_MAXPRIMERSATSITE="10"
ARG_MAXTARGETLENGTH="500"
ARG_OBJFNWEIGHTS="0.50 0.25"
ARG_BESTNTARGETS="30"
ARG_PREDICTIVE_MODELS="${PREDICTIVE_MODELS_PATH}/classify/model-51373185 ${PREDICTIVE_MODELS_PATH}/regress/model-f8b6fd5d"


# Make tmp directory for memoizing alignments and stats
mkdir -p $PREP_MEMOIZE_DIR

function run_for_taxid() {
    # Set information on taxonomy, from arguments
    taxid="$1"
    segment="$2"
    refaccs="$3"

    echo "Running for taxid $taxid (segment: $segment)" > /dev/tty

    # Make an output directory
    outdir="tax-${taxid}_${segment}"
    mkdir -p $outdir

    # Fetch accessions and create table of them
    conda activate adapt
    if [ ! -f $outdir/accessions.tsv ]; then
        python ../scripts/find_year_for_accessions.py $taxid $segment | awk -v taxid="$taxid" -v segment="$segment" '{print taxid"\t"segment"\t"$1}' | sort | uniq > $outdir/accessions.tsv
    fi

    # Randomly sample a number of sequences equal to the number
    # of accessions (bootstrapping over the input)
    SAMPLE_SIZE=$(cat $outdir/accessions.tsv | wc -l)

    # Produce NUM_DESIGNS designs, in parallel
    mkdir -p $outdir/designs/resampled
    mkdir -p $outdir/designs/non-resampled

    function design() {
        # Echo to /dev/tty so that it's printed without waiting for the function to end
        echo "Starting design for $taxid (segment: $segment), $1 of $NUM_DESIGNS" > /dev/tty

        # Sleep 0-60 seconds, so there are not too many NCBI requests at once
        sleep $((RANDOM % 60))

        if [ ! -f $outdir/designs/resampled/design-${1}.tsv.0 ]; then
            design.py complete-targets auto-from-args $taxid "$segment" $refaccs $outdir/designs/resampled/design-${1}.tsv --obj maximize-activity --soft-guide-constraint $ARG_SOFTGUIDECONSTRAINT --hard-guide-constraint $ARG_HARDGUIDECONSTRAINT --penalty-strength $ARG_PENALTYSTRENGTH --maximization-algorithm $ARG_MAXIMIZATIONALGORITHM -gl $ARG_GL -pl $ARG_PL -pm $ARG_PM -pp $ARG_PP --primer-gc-content-bounds $ARG_PRIMER_GC_LO $ARG_PRIMER_GC_HI --max-primers-at-site $ARG_MAXPRIMERSATSITE --max-target-length $ARG_MAXTARGETLENGTH --obj-fn-weights $ARG_OBJFNWEIGHTS --best-n-targets $ARG_BESTNTARGETS --predict-activity-model-path $ARG_PREDICTIVE_MODELS --mafft-path $MAFFT_PATH --prep-memoize-dir $PREP_MEMOIZE_DIR --ncbi-api-key $NCBI_API_KEY --sample-seqs $SAMPLE_SIZE --cluster-threshold $CLUSTER_THRESHOLD --use-accessions $outdir/accessions.tsv --verbose &> $outdir/designs/resampled/design-${1}.out
        else
            echo "Resampled design for $taxid (segment: $segment), $1 of $NUM_DESIGNS, already exists; skipping" > /dev/tty
        fi
        if [ ! -f $outdir/designs/non-resampled/design-${1}.tsv.0 ]; then
            design.py complete-targets auto-from-args $taxid "$segment" $refaccs $outdir/designs/non-resampled/design-${1}.tsv --obj maximize-activity --soft-guide-constraint $ARG_SOFTGUIDECONSTRAINT --hard-guide-constraint $ARG_HARDGUIDECONSTRAINT --penalty-strength $ARG_PENALTYSTRENGTH --maximization-algorithm $ARG_MAXIMIZATIONALGORITHM -gl $ARG_GL -pl $ARG_PL -pm $ARG_PM -pp $ARG_PP --primer-gc-content-bounds $ARG_PRIMER_GC_LO $ARG_PRIMER_GC_HI --max-primers-at-site $ARG_MAXPRIMERSATSITE --max-target-length $ARG_MAXTARGETLENGTH --obj-fn-weights $ARG_OBJFNWEIGHTS --best-n-targets $ARG_BESTNTARGETS --predict-activity-model-path $ARG_PREDICTIVE_MODELS --mafft-path $MAFFT_PATH --prep-memoize-dir $PREP_MEMOIZE_DIR --ncbi-api-key $NCBI_API_KEY --cluster-threshold $CLUSTER_THRESHOLD --use-accessions $outdir/accessions.tsv --verbose &> $outdir/designs/non-resampled/design-${1}.out
        else
            echo "Non-resampled design for $taxid (segment: $segment), $1 of $NUM_DESIGNS, already exists; skipping" > /dev/tty
        fi
        echo "Completed design for $taxid (segment: $segment), $1 of $NUM_DESIGNS" > /dev/tty
    }

    # Export variables, so they can be read by parallel
    export taxid
    export segment
    export refaccs
    export outdir
    export NUM_DESIGNS
    export PREP_MEMOIZE_DIR
    export MAFFT_PATH
    export SAMPLE_SIZE
    export CLUSTER_THRESHOLD
    export NCBI_API_KEY
    export ARG_GL
    export ARG_PL
    export ARG_PM
    export ARG_PP
    export ARG_PRIMER_GC_LO
    export ARG_PRIMER_GC_HI
    export ARG_SOFTGUIDECONSTRAINT
    export ARG_HARDGUIDECONSTRAINT
    export ARG_PENALTYSTRENGTH
    export ARG_MAXIMIZATIONALGORITHM
    export ARG_MAXPRIMERSATSITE
    export ARG_MAXTARGETLENGTH
    export ARG_OBJFNWEIGHTS
    export ARG_BESTNTARGETS
    export ARG_PREDICTIVE_MODELS
    export -f design

    # Run parallel
    parallel --jobs $NJOBS --no-notice design ::: $(seq 1 $NUM_DESIGNS)

    echo "Done running for taxid $taxid (segment: $segment)" > /dev/tty
}


# Run for Zika virus
run_for_taxid "64320" "None" "NC_035889,NC_012532"

# Run for Lassa virus, S segment
run_for_taxid "11620" "S" "KM821998,GU481072,KM821773"

# Run for Lassa virus, L segment
run_for_taxid "11620" "L" "U73034"

# Run for Ebola virus (Zaire)
#run_for_taxid "186538" "None" "NC_002549"

# Run for Nipah virus
#run_for_taxid "121791" "None" "NC_002728"

# Run for HIV-1
run_for_taxid "11676" "None" "NC_001802"

# Run for HCV
#run_for_taxid "11103" "None" "NC_004102,NC_030791,NC_009827,NC_009826,NC_009825,NC_038882,NC_009824,NC_009823"

# Run for IAV segment 2
run_for_taxid "11320" "2" "NC_026435,NC_002021,NC_007375,NC_026423,NC_007372"

# Run for Rhinovirus A
run_for_taxid "147711" "None" "NC_038311,NC_001617,NC_038311"

# Run for Enterovirus A
run_for_taxid "138948" "None" "NC_038306,NC_001612,NC_038306"
