#!/bin/bash

# Make designs to use for comparing performance against baseline (naive) designs.
#
# To make a more clear comparison to naive designs (which design one guide per
# window), this only uses adapt with the sliding window approach.
#
# This decides that a guide hits (detects) a target if it is within 1 mismatch
# (since ARG_GM is 1). Note that this *does* tolerate G-U pairs -- G-U pairs are
# not counted as mismatches. If we did not want to count G-U pairs, we would
# use --do-not-allow-gu-pairs.
#
# Author: Hayden Metsky <hayden@mit.edu>


# Load environment, and variables, for ADAPT
source ~/misc-repos/adapt-designs/scripts/run-adapt/custom-env/load_custom_env.sh

# Set variables for measuring uncertainty
NUM_DESIGNS=5

# Set variables for running
NJOBS=8

# Set variables for design
CLUSTER_THRESHOLD=0.4   # Use high value to make it more likely to obtain only one cluster or a large one
ARG_WINDOWSIZE="200"
ARG_GL="30"
ARG_GM="1"

# For minimize-guides
ARG_GP="0.99"

# For maximize-activity; note that this will use --use-simple-binary-activity-prediction
# We will vary the hard guide constraint; set the penalty strength to 0 to ignore
# the soft guide constraint
ARG_MAXIMIZATIONALGORITHM="random-greedy"
ARG_PENALTYSTRENGTH="0"

# Write commands to a file
commands_fn=$(mktemp)
echo -n "" > $commands_fn

function run_for_taxid() {
    # Set information on taxonomy, from arguments
    taxid="$1"
    segment="$2"
    refaccs="$3"

    echo "Creating commands for taxid $taxid (segment: $segment)" > /dev/tty

    # Make an output directory
    outdir="tax-${taxid}_${segment}"
    mkdir -p $outdir
    mkdir -p $outdir/input-alns
    mkdir -p $outdir/designs

    conda activate adapt

    # Fetch accessions and create table of them
    if [ ! -f $outdir/accessions.tsv ]; then
        python ../scripts/find_year_for_accessions.py $taxid $segment | awk -v taxid="$taxid" -v segment="$segment" '{print taxid"\t"segment"\t"$1}' | sort | uniq > $outdir/accessions.tsv
    fi

    # Run design.py to create an alignment from all the sequences; continue with
    # the sliding-window design approach (using a window size equal to guide size
    # so it's fast), but ignore the output
    # Only do so if an alignment does not already exist
    if [ ! -f $outdir/input-alns/all-accessions.fasta ]; then
        design.py sliding-window auto-from-args $taxid $segment $refaccs /tmp/design-for-aln.tsv --window-size $ARG_GL -gl $ARG_GL -gm $ARG_GM -gp $ARG_GP --mafft-path $MAFFT_PATH --prep-memoize-dir $PREP_MEMOIZE_DIR --cluster-threshold $CLUSTER_THRESHOLD --use-accessions $outdir/accessions.tsv --write-input-aln $outdir/input-alns/all-accessions.fasta --verbose &> $outdir/input-alns/all-accessions.out
        rm /tmp/design-for-aln.tsv.*
        # An output alignment will have a `.0` suffix because it corresponds to the
        # first cluster; rename it to remove the suffix and delete all other clusters
        # (i.e., only use the first, which is the largest)
        mv $outdir/input-alns/all-accessions.fasta.0 $outdir/input-alns/all-accessions.fasta
        rm -f $outdir/input-alns/all-accessions.fasta.*
    fi

    # Randomly sample a number of sequences equal to the number
    # of sequences in the FASTA (bootstrapping over the input)
    sample_size=$(cat $outdir/input-alns/all-accessions.fasta | grep '>' | wc -l)

    # Randomly sample the alignment to create a fasta file as input for
    # each design; by sampling the same alignment, all designs are
    # produced from the same coordinate space and so window positions
    # across designs can be matched
    # (If we were to use --sample-seqs with design.py, it would create a
    # separate alignment for each design where the coordinates could
    # differ across those alignments)
    for i in $(seq 1 $NUM_DESIGNS); do
        python ../scripts/randomly_sample_fasta.py $outdir/input-alns/all-accessions.fasta $sample_size $outdir/input-alns/design-${i}.fasta
    done

    # Produce a design.py command for each design, using the alignment
    # produced above, with the minimize-guides objective
    for i in $(seq 1 $NUM_DESIGNS); do
        echo "design.py sliding-window fasta $outdir/input-alns/design-${i}.fasta -o $outdir/designs/design-${i}.real-design.minimize-guides.tsv --window-size $ARG_WINDOWSIZE -gl $ARG_GL --obj minimize-guides -gm $ARG_GM -gp $ARG_GP --verbose &> $outdir/designs/design-${i}.real-design.minimize-guides.out" >> $commands_fn
    done

    # Do the same with the maximize-activity objective and different hard guide
    # constraints
    for i in $(seq 1 $NUM_DESIGNS); do
        for hgc in 1 2 3 4 5; do
            echo "design.py sliding-window fasta $outdir/input-alns/design-${i}.fasta -o $outdir/designs/design-${i}.real-design.maximize-activity.hgc-${hgc}.tsv --window-size $ARG_WINDOWSIZE -gl $ARG_GL --obj maximize-activity -gm $ARG_GM --hard-guide-constraint $hgc --use-simple-binary-activity-prediction --maximization-algorithm $ARG_MAXIMIZATIONALGORITHM --penalty-strength $ARG_PENALTYSTRENGTH --verbose &> $outdir/designs/design-${i}.real-design.maximize-activity.hgc-${hgc}.out" >> $commands_fn
        done
    done

    # Produce a design_naively.py command for each design, using the
    # alignment produced above
    for i in $(seq 1 $NUM_DESIGNS); do
        echo "design_naively.py $outdir/input-alns/design-${i}.fasta $outdir/designs/design-${i}.naive-design.tsv --window-size $ARG_WINDOWSIZE -gl $ARG_GL -gm $ARG_GM --verbose &> $outdir/designs/design-${i}.naive-design.out" >> $commands_fn
    done
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
run_for_taxid "11103" "None" "NC_004102,NC_030791,NC_009827,NC_009826,NC_009825,NC_038882,NC_009824,NC_009823"

# Run for IAV segment 2
run_for_taxid "11320" "2" "NC_026435,NC_002021,NC_007375,NC_026423,NC_007372"

# Run for human coronavirus 229E
#run_for_taxid "11137" "None" "NC_028752"

# Run for Rhinovirus A, B, C
run_for_taxid "147711" "None" "NC_038311,NC_001617,NC_038311"
run_for_taxid "147712" "None" "NC_038312"
run_for_taxid "463676" "None" "NC_038878"

# Run for Enterovirus A, B, C, D
run_for_taxid "138948" "None" "NC_038306,NC_001612,NC_038306"
run_for_taxid "138949" "None" "NC_001472,NC_038307,NC_038307,NC_001472"
run_for_taxid "138950" "None" "NC_002058"
run_for_taxid "138951" "None" "NC_038308,NC_001430"

# Run parallel on the commands
echo "Running all commands.."
parallel --jobs $NJOBS --no-notice < $commands_fn
rm $commands_fn
