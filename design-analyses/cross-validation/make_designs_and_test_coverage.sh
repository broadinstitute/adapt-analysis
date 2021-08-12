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
NJOBS=20
CLUSTER_THRESHOLD=1.0   # Use high value to obtain a single cluster
ARG_GL="28"
ARG_PL="30"
ARG_PM_DESIGN="3"
ARG_PM_ANALYZE="3"
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
ARG_PREDICTIVE_MODELS="${PREDICTIVE_MODELS_PATH}/classify/cas13a/v1_0 ${PREDICTIVE_MODELS_PATH}/regress/cas13a/v1_0"

# Set some more relaxed variables (they are relaxed for the guide activity
# objective - i.e., relaxed constraints to achieve higher activity - but are
# actually more strict for primers in order to achieve higher coverage)
ARG_PM_DESIGN_RELAXED="3"   # same as standard
ARG_PM_ANALYZE_RELAXED="4"
ARG_PP_RELAXED="0.995"
ARG_PRIMER_GC_LO_RELAXED="0.20"
ARG_PRIMER_GC_HI_RELAXED="0.80"
ARG_SOFTGUIDECONSTRAINT_RELAXED="3"
ARG_HARDGUIDECONSTRAINT_RELAXED="10"
ARG_PENALTYSTRENGTH_RELAXED="0.05"
ARG_MAXPRIMERSATSITE_RELAXED="15"
ARG_MAXTARGETLENGTH_RELAXED="1000"
ARG_OBJFNWEIGHTS_RELAXED="0.30 0.05"

# Make tmp directory for memoizing alignments and stats
mkdir -p $PREP_MEMOIZE_DIR

# Set tmp directory
if [ -d "/ebs/tmpfs/tmp" ]; then
    export TMPDIR="/ebs/tmpfs/tmp"
else
    export TMPDIR="/tmp"
fi


function run_for_taxid() {
    run_for_taxid_and_constrainttype $1 $2 $3 "standard"
    run_for_taxid_and_constrainttype $1 $2 $3 "relaxed"
}

function run_for_taxid_and_constrainttype() {
    # Set information on taxonomy, from arguments
    taxid="$1"
    segment="$2"
    refaccs="$3"

    echo "Running for taxid $taxid (segment: $segment)" > /dev/tty

    # ADAPT's --ref-accs wants the reference accessions to be space-separated
    refaccs=$(echo "$refaccs" | tr ',' ' ')

    # Set the constraint type, from arguments ('standard' or 'relaxed')
    constrainttype="$4"
    if [ "$constrainttype" == "standard" ]; then
        echo "With standard constraints on guides" > /dev/tty
        sgc="$ARG_SOFTGUIDECONSTRAINT"
        hgc="$ARG_HARDGUIDECONSTRAINT"
        ps="$ARG_PENALTYSTRENGTH"
        pmdesign="$ARG_PM_DESIGN"
        pmanalyze="$ARG_PM_ANALYZE"
        pp="$ARG_PP"
        gclo="$ARG_PRIMER_GC_LO"
        gchi="$ARG_PRIMER_GC_HI"
        mps="$ARG_MAXPRIMERSATSITE"
        mtl="$ARG_MAXTARGETLENGTH"
        objfnweights="$ARG_OBJFNWEIGHTS"
        analyzepredictthresargs=""  # use the default for the model
    elif [ "$constrainttype" == "relaxed" ]; then
        echo "With relaxed constraints on guides" > /dev/tty
        sgc="$ARG_SOFTGUIDECONSTRAINT_RELAXED"
        hgc="$ARG_HARDGUIDECONSTRAINT_RELAXED"
        ps="$ARG_PENALTYSTRENGTH_RELAXED"
        pmdesign="$ARG_PM_DESIGN_RELAXED"
        pmanalyze="$ARG_PM_ANALYZE_RELAXED"
        pp="$ARG_PP_RELAXED"
        gclo="$ARG_PRIMER_GC_LO_RELAXED"
        gchi="$ARG_PRIMER_GC_HI_RELAXED"
        mps="$ARG_MAXPRIMERSATSITE_RELAXED"
        mtl="$ARG_MAXTARGETLENGTH_RELAXED"
        objfnweights="$ARG_OBJFNWEIGHTS_RELAXED"
        analyzepredictthresargs="--predict-activity-thres 0.3 2.7198637"   # arbitrary for classification (0.3; lower precision but higher sensitivity than default) and default for regression (-1.28 + 4)
    else
        echo "Unknown constraint type '$constrainttype'"
        exit 1
    fi

    # Make an output directory
    outdir="tax-${taxid}_${segment}"
    mkdir -p $outdir

    mkdir -p $outdir/designs
    mkdir -p $outdir/designs/accessions
    mkdir -p $outdir/designs/$constrainttype/designs
    mkdir -p $outdir/designs/$constrainttype/coverages

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

        if [ ! -f $outdir/designs/$constrainttype/designs/design-${i}.tsv.0 ]; then
            # Produce a design command
            echo -n "design.py complete-targets auto-from-args $taxid $segment $outdir/designs/$constrainttype/designs/design-${i}.tsv --obj maximize-activity --soft-guide-constraint $sgc --hard-guide-constraint $hgc --penalty-strength $ps --maximization-algorithm $ARG_MAXIMIZATIONALGORITHM -gl $ARG_GL -pl $ARG_PL -pm $pmdesign -pp $pp --primer-gc-content-bounds $gclo $gchi --max-primers-at-site $mps --max-target-length $mtl --obj-fn-weights $objfnweights --best-n-targets $ARG_BESTNTARGETS --predict-activity-model-path $ARG_PREDICTIVE_MODELS --mafft-path $MAFFT_PATH --prep-memoize-dir $PREP_MEMOIZE_DIR --ncbi-api-key $NCBI_API_KEY --cluster-threshold $CLUSTER_THRESHOLD --use-accessions $outdir/designs/accessions/design-${i}.design-set.tsv --ref-accs $refaccs --verbose &> $outdir/designs/$constrainttype/designs/design-${i}.out" >> $commands_fn
            echo -n "; " >> $commands_fn
        fi

        if [ ! -f $outdir/designs/$constrainttype/coverages/design-${i}.coverage-against-test.active.txt ]; then
            # Produce analyze_coverage command, with decisions made based on whether pairs are active
            cat $outdir/designs/accessions/design-${i}.test-set.tsv | awk '{print $3}' > $outdir/designs/accessions/design-${i}.test-set.acc-only.txt
            echo -n "analyze_coverage.py $outdir/designs/$constrainttype/designs/design-${i}.tsv.0 $outdir/designs/accessions/design-${i}.test-set.acc-only.txt --predict-activity-model-path $ARG_PREDICTIVE_MODELS $analyzepredictthresargs -pm $pmanalyze --use-accessions --write-frac-bound $outdir/designs/$constrainttype/coverages/design-${i}.coverage-against-test.active.txt --verbose &> $outdir/designs/$constrainttype/coverages/design-${i}.coverage-against-test.active.out" >> $commands_fn
            echo -n "; " >> $commands_fn
        fi

        if [ ! -f $outdir/designs/$constrainttype/coverages/design-${i}.coverage-against-test.highly-active.txt ]; then
            # Produce analyze_coverage command, with decisions made based on whether pairs are *highly* active
            cat $outdir/designs/accessions/design-${i}.test-set.tsv | awk '{print $3}' > $outdir/designs/accessions/design-${i}.test-set.acc-only.txt
            echo -n "analyze_coverage.py $outdir/designs/$constrainttype/designs/design-${i}.tsv.0 $outdir/designs/accessions/design-${i}.test-set.acc-only.txt --predict-activity-model-path $ARG_PREDICTIVE_MODELS $analyzepredictthresargs --predict-activity-require-highly-active -pm $pmanalyze --use-accessions --write-frac-bound $outdir/designs/$constrainttype/coverages/design-${i}.coverage-against-test.highly-active.txt --verbose &> $outdir/designs/$constrainttype/coverages/design-${i}.coverage-against-test.highly-active.out" >> $commands_fn
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
##run_for_taxid "11103" "None" "NC_004102,NC_030791,NC_009827,NC_009826,NC_009825,NC_038882,NC_009824,NC_009823"

# Run for IAV segment 2
##run_for_taxid "11320" "2" "NC_026435,NC_002021,NC_007375,NC_026423,NC_007372"

# Run for Rhinovirus A
run_for_taxid "147711" "None" "NC_038311,NC_001617,NC_038311"

# Run for Enterovirus A
run_for_taxid "138948" "None" "NC_038306,NC_001612,NC_038306"
