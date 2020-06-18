#!/bin/bash

# Make designs on a subset of the input data (the 'design set') and measure the
# coverage these changes achieve against the remaining data (the 'test set').
# We can think of this like performing Monte Carlo cross-validation (or
# repeated randomly subsampling validation) of the design method.
#
# Author: Hayden Metsky <hayden@mit.edu>


# Load environment, and variables, for ADAPT
source ~/misc-repos/adapt-designs/scripts/run-adapt/custom-env/load_custom_env.sh

# Set variables for the cross-validation
NUM_DESIGNS=20
DESIGN_SET_FRACTION=0.8 # Use 80% of accessions as input for design; test against remaining 20%

# Set variables for design
NJOBS=10
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

# Set tmp directory
if [ -d "/ebs/tmpfs/tmp" ]; then
    export TMPDIR="/ebs/tmpfs/tmp"
else
    export TMPDIR="/tmp"
fi


function run_for_taxid() {
    # Set information on taxonomy, from arguments
    taxid="$1"
    segment="$2"
    refaccs="$3"

    echo "Running for taxid $taxid (segment: $segment)" > /dev/tty

    # Make an output directory
    outdir="tax-${taxid}_${segment}"
    mkdir -p $outdir

    mkdir -p $outdir/designs
    mkdir -p $outdir/designs/accessions
    mkdir -p $outdir/designs/designs
    mkdir -p $outdir/designs/coverages

    conda activate adapt

    # Fetch accessions and create table of them
    conda activate adapt
    if [ ! -f $outdir/accessions.tsv ]; then
        python ../scripts/find_year_for_accessions.py $taxid $segment | awk -v taxid="$taxid" -v segment="$segment" '{print taxid"\t"segment"\t"$1}' | sort | uniq > $outdir/accessions.tsv
    fi

    # Determine number of accessions in the design set (rounded to nearest integer)
    total_num_acc=$(cat $outdir/accessions.tsv | wc -l)
    design_set_num_acc=$(echo "$total_num_acc" | awk -v frac="$DESIGN_SET_FRACTION" '{print int($1*frac + 0.5)}')
    test_set_num_acc=$(echo "$total_num_acc" | awk -v designsetnum="$design_set_num_acc" '{print $1 - designsetnum}')

    # Write commands to a file
    commands_fn=$(mktemp)
    echo -n "" > $commands_fn

    # Make commands, with $NUM_DESIGNS splits of the data (each with randomly split design/test data sets)
    for i in $(seq 1 $NUM_DESIGNS); do
        if [ ! -f $outdir/designs/accessions/design-${i}.test-set.tsv ]; then
            # Split all the accessions into a set for design and a set for testing
            sort -R $outdir/accessions.tsv > $outdir/accessions.shuffled.tsv
            head -n $design_set_num_acc $outdir/accessions.shuffled.tsv > $outdir/designs/accessions/design-${i}.design-set.tsv
            tail -n $test_set_num_acc $outdir/accessions.shuffled.tsv > $outdir/designs/accessions/design-${i}.test-set.tsv
            rm $outdir/accessions.shuffled.tsv
        fi

        # Produce a design.py and analyze_coverage.py command

        if [ ! -f $outdir/designs/designs/design-${i}.tsv.0 ]; then
            # Produce a design command
            echo -n "design.py complete-targets auto-from-args $taxid $segment $refaccs $outdir/designs/designs/design-${i}.tsv --obj maximize-activity --soft-guide-constraint $ARG_SOFTGUIDECONSTRAINT --hard-guide-constraint $ARG_HARDGUIDECONSTRAINT --penalty-strength $ARG_PENALTYSTRENGTH --maximization-algorithm $ARG_MAXIMIZATIONALGORITHM -gl $ARG_GL -pl $ARG_PL -pm $ARG_PM -pp $ARG_PP --primer-gc-content-bounds $ARG_PRIMER_GC_LO $ARG_PRIMER_GC_HI --max-primers-at-site $ARG_MAXPRIMERSATSITE --max-target-length $ARG_MAXTARGETLENGTH --obj-fn-weights $ARG_OBJFNWEIGHTS --best-n-targets $ARG_BESTNTARGETS --predict-activity-model-path $ARG_PREDICTIVE_MODELS --mafft-path $MAFFT_PATH --prep-memoize-dir $PREP_MEMOIZE_DIR --ncbi-api-key $NCBI_API_KEY --cluster-threshold $CLUSTER_THRESHOLD --use-accessions $outdir/designs/accessions/design-${i}.design-set.tsv --verbose &> $outdir/designs/designs/design-${i}.out" >> $commands_fn
            echo -n "; " >> $commands_fn
        fi

        if [ ! -f $outdir/designs/coverages/design-${i}.coverage-against-test.active.txt ]; then
            # Produce analyze_coverage command, with decisions made based on whether pairs are active
            cat $outdir/designs/accessions/design-${i}.test-set.tsv | awk '{print $3}' > $outdir/designs/accessions/design-${i}.test-set.acc-only.txt
            echo -n "analyze_coverage.py $outdir/designs/designs/design-${i}.tsv.0 $outdir/designs/accessions/design-${i}.test-set.acc-only.txt --predict-activity-model-path $ARG_PREDICTIVE_MODELS -pm $ARG_PM --use-accessions --write-frac-bound $outdir/designs/coverages/design-${i}.coverage-against-test.active.txt --verbose &> $outdir/designs/coverages/design-${i}.coverage-against-test.active.out" >> $commands_fn
            echo "" >> $commands_fn
        fi

        if [ ! -f $outdir/designs/coverages/design-${i}.coverage-against-test.highly-active.txt ]; then
            # Produce analyze_coverage command, with decisions made based on whether pairs are *highly* active
            cat $outdir/designs/accessions/design-${i}.test-set.tsv | awk '{print $3}' > $outdir/designs/accessions/design-${i}.test-set.acc-only.txt
            echo -n "analyze_coverage.py $outdir/designs/designs/design-${i}.tsv.0 $outdir/designs/accessions/design-${i}.test-set.acc-only.txt --predict-activity-model-path $ARG_PREDICTIVE_MODELS --predict-activity-require-highly-active -pm $ARG_PM --use-accessions --write-frac-bound $outdir/designs/coverages/design-${i}.coverage-against-test.highly-active.txt --verbose &> $outdir/designs/coverages/design-${i}.coverage-against-test.highly-active.out" >> $commands_fn
            echo "" >> $commands_fn
        fi
    done

    # Run parallel
    parallel --jobs $NJOBS --no-notice < $commands_fn

    rm $commands_fn
}


# Run for Zika virus
run_for_taxid "64320" "None" "NC_035889,NC_012532"

# Run for Lassa virus, S segment
run_for_taxid "11620" "S" "KM821998,GU481072,KM821773"

# Run for Lassa virus, L segment
run_for_taxid "11620" "L" "U73034"

# Run for Ebola virus (Zaire)
run_for_taxid "186538" "None" "NC_002549"

# Run for Nipah virus
run_for_taxid "121791" "None" "NC_002728"

# Run for HIV-1
run_for_taxid "11676" "None" "NC_001802"

# Run for HCV
#run_for_taxid "11103" "None" "NC_004102,NC_030791,NC_009827,NC_009826,NC_009825,NC_038882,NC_009824,NC_009823"

# Run for IAV segment 2
#run_for_taxid "11320" "2" "NC_026435,NC_002021,NC_007375,NC_026423,NC_007372"

# Run for Rhinovirus A
run_for_taxid "147711" "None" "NC_038311,NC_001617,NC_038311"

# Run for Enterovirus A
run_for_taxid "138948" "None" "NC_038306,NC_001612,NC_038306"
