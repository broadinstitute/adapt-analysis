#!/bin/bash

# Make designs to use for evaluating dispersion.
#
# Author: Hayden Metsky <hayden@mit.edu>


# Allow activating conda environments
source ~/anaconda3/etc/profile.d/conda.sh

# Set variables for measuring dispersion
NUM_DESIGNS=100

# Set variables for design
NJOBS=8
PREP_MEMOIZE_DIR="/ebs/dgd-analysis/prep-memoize-dir"
MAFFT_PATH="/home/hayden/viral-ngs/viral-ngs-etc/conda-env/bin/mafft"
CLUSTER_THRESHOLD=1.0   # Use high value to obtain a single cluster
ARG_GL="28"
ARG_GM="1"
ARG_GP="0.99"
ARG_PL="30"
ARG_PM="3"
ARG_PP="0.99"
ARG_MAXPRIMERSATSITE="10"
ARG_MAXTARGETLENGTH="1000"
ARG_COSTFNWEIGHTS="0.6667 0.2222 0.1111"
ARG_BESTNTARGETS="30"


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
    conda activate data-analysis
    if [ ! -f $outdir/accessions.tsv ]; then
        python ../scripts/find_year_for_accessions.py $taxid $segment | awk -v taxid="$taxid" -v segment="$segment" '{print taxid"\t"segment"\t"$1}' | sort | uniq > $outdir/accessions.tsv
    fi

    # Randomly sample a number of sequences equal to the number
    # of accessions (bootstrapping over the input)
    SAMPLE_SIZE=$(cat $outdir/accessions.tsv | wc -l)

    # Produce NUM_DESIGNS designs, in parallel
    conda activate dgd
    mkdir -p $outdir/designs/resampled
    mkdir -p $outdir/designs/non-resampled

    function design() {
        # Echo to /dev/tty so that it's printed without waiting for the function to end
        echo "Starting design for $taxid (segment: $segment), $1 of $NUM_DESIGNS" > /dev/tty

        # Sleep 0-60 seconds, so there are not too many NCBI requests at once
        sleep $((RANDOM % 60))

        design.py complete-targets auto-from-args $taxid $segment $refaccs $outdir/designs/resampled/design-${1}.tsv -gl $ARG_GL -gm $ARG_GM -gp $ARG_GP -pl $ARG_PL -pm $ARG_PM -pp $ARG_PP --max-primers-at-site $ARG_MAXPRIMERSATSITE --max-target-length $ARG_MAXTARGETLENGTH --cost-fn-weights $ARG_COSTFNWEIGHTS --best-n-targets $ARG_BESTNTARGETS --mafft-path $MAFFT_PATH --prep-memoize-dir $PREP_MEMOIZE_DIR --sample-seqs $SAMPLE_SIZE --cluster-threshold $CLUSTER_THRESHOLD --use-accessions $outdir/accessions.tsv --verbose &> $outdir/designs/resampled/design-${1}.out
        design.py complete-targets auto-from-args $taxid $segment $refaccs $outdir/designs/non-resampled/design-${1}.tsv -gl $ARG_GL -gm $ARG_GM -gp $ARG_GP -pl $ARG_PL -pm $ARG_PM -pp $ARG_PP --max-primers-at-site $ARG_MAXPRIMERSATSITE --max-target-length $ARG_MAXTARGETLENGTH --cost-fn-weights $ARG_COSTFNWEIGHTS --best-n-targets $ARG_BESTNTARGETS --mafft-path $MAFFT_PATH --prep-memoize-dir $PREP_MEMOIZE_DIR --cluster-threshold $CLUSTER_THRESHOLD --use-accessions $outdir/accessions.tsv --verbose &> $outdir/designs/non-resampled/design-${1}.out
        echo "Completed design for $taxid (segment: $segment), $1 of $NUM_DESIGNS" > /dev/tty
    }

    # Export variables, so they can be read by parallel
    export taxid
    export segment
    export refaccs
    export outdir
    export NUM_DESIGNS
    export MAFFT_PATH
    export SAMPLE_SIZE
    export CLUSTER_THRESHOLD
    export ARG_GL
    export ARG_GM
    export ARG_GP
    export ARG_PL
    export ARG_PM
    export ARG_PP
    export ARG_MAXPRIMERSATSITE
    export ARG_MAXTARGETLENGTH
    export ARG_COSTFNWEIGHTS
    export ARG_BESTNTARGETS
    export -f design

    # Run parallel
    parallel --jobs $NJOBS --no-notice design ::: $(seq 1 $NUM_DESIGNS)

    echo "Done running for taxid $taxid (segment: $segment)" > /dev/tty
}


# Run for Zika virus
run_for_taxid "64320" "None" "NC_035889,NC_012532"

# Run for Lassa virus, S segment
run_for_taxid "11620" "S" "KM821998,GU481072,KM821773"
