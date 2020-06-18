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


# Allow activating conda environments
source ~/anaconda3/etc/profile.d/conda.sh

# Set variables for measuring uncertainty
NUM_DESIGNS=5

# Set variables for running
NJOBS=8

# Set variables for design
CLUSTER_THRESHOLD=0.4   # Use high value to make it more likely to obtain only one cluster or a large one
ARG_WINDOWSIZE="200"
ARG_GL="30"
ARG_GM="1"

# For maximize-activity; note that this will use --use-simple-binary-activity-prediction
# We will vary the hard guide constraint; set the penalty strength to 0 to ignore
# the soft guide constraint
ARG_MAXIMIZATIONALGORITHM="random-greedy"
ARG_PENALTYSTRENGTH="0"

# Set tmp directory
if [ -d "/ebs/tmpfs/tmp" ]; then
    export TMPDIR="/ebs/tmpfs/tmp"
else
    export TMPDIR="/tmp"
fi

# Write commands to a file
commands_fn=$(mktemp)
echo -n "" > $commands_fn

function run_for_taxid() {
    # Set information on taxonomy, from arguments
    taxid="$1"
    segment="$2"
    refaccs="$3"
    alnref="$4"

    echo "Creating commands for taxid $taxid (segment: $segment)" > /dev/tty

    # Make an output directory
    outdir="tax-${taxid}_${segment}"
    mkdir -p $outdir
    mkdir -p $outdir/input-alns
    mkdir -p $outdir/designs

    # Load adapt environment
    conda activate adapt

    # Fetch accessions and create table of them
    if [ ! -f $outdir/accessions.tsv ]; then
        python ../scripts/find_year_for_accessions.py $taxid $segment | awk -v taxid="$taxid" -v segment="$segment" '{print taxid"\t"segment"\t"$1}' | sort | uniq > $outdir/accessions.tsv
    fi

    # Create an alignment against a reference using mafft --addfragments --keeplength; this
    # creates a more clean MSA without many gaps, but does remove sequence from genomes
    # if they are insertions relative to the reference
    # Only do so if an alignment does not already exist
    if [ ! -f $outdir/input-alns/all-accessions.fasta ]; then
        # Split accessions into file with only $alnref and file without $alnref
        accessions_alnref=$(mktemp)
        echo "$alnref" > $accessions_alnref
        accessions_other=$(mktemp)
        cat $outdir/accessions.tsv | awk -v alnref="$alnref" '$3 != alnref {print $3}' > $accessions_other

        # Download FASTA for the above accession files
        fasta_alnref=$(mktemp)
        fasta_other=$(mktemp)
        python ../scripts/download_fasta_for_accessions.py $accessions_alnref > $fasta_alnref
        python ../scripts/download_fasta_for_accessions.py $accessions_other > $fasta_other

        # Load environment with mafft
        conda activate misc-bioinformatic-tools

        # Align with mafft
        mafft --preservecase --auto --keeplength --addfragments $fasta_other $fasta_alnref > $outdir/input-alns/all-accessions.fasta 2> $outdir/input-alns/all-accessions.out

        # Reload adapt
        conda activate adapt

        rm $accessions_alnref
        rm $accessions_other
        rm $fasta_alnref
        rm $fasta_other
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
        if [ ! -f $outdir/input-alns/design-${i}.fasta ]; then
            python ../scripts/randomly_sample_fasta.py $outdir/input-alns/all-accessions.fasta $sample_size $outdir/input-alns/design-${i}.fasta
        fi
    done

    # Produce a design.py command for each design, using the alignment
    # produced above, with the minimize-guides objective
    for i in $(seq 1 $NUM_DESIGNS); do
        for gp in 0.9 0.95 0.99; do
            if [ ! -f $outdir/designs/design-${i}.real-design.minimize-guides.gp-${gp}.tsv ]; then
                echo "design.py sliding-window fasta $outdir/input-alns/design-${i}.fasta -o $outdir/designs/design-${i}.real-design.minimize-guides.gp-${gp}.tsv --window-size $ARG_WINDOWSIZE -gl $ARG_GL --obj minimize-guides -gm $ARG_GM -gp $gp --verbose &> $outdir/designs/design-${i}.real-design.minimize-guides.gp-${gp}.out" >> $commands_fn
            fi
        done
    done

    # Do the same with the maximize-activity objective and different hard guide
    # constraints
    for i in $(seq 1 $NUM_DESIGNS); do
        for hgc in 1 2 3 4 5; do
            if [ ! -f $outdir/designs/design-${i}.real-design.maximize-activity.hgc-${hgc}.tsv ]; then
                echo "design.py sliding-window fasta $outdir/input-alns/design-${i}.fasta -o $outdir/designs/design-${i}.real-design.maximize-activity.hgc-${hgc}.tsv --window-size $ARG_WINDOWSIZE -gl $ARG_GL --obj maximize-activity -gm $ARG_GM --hard-guide-constraint $hgc --use-simple-binary-activity-prediction --maximization-algorithm $ARG_MAXIMIZATIONALGORITHM --penalty-strength $ARG_PENALTYSTRENGTH --verbose &> $outdir/designs/design-${i}.real-design.maximize-activity.hgc-${hgc}.out" >> $commands_fn
            fi
        done
    done

    # Produce a design_naively.py command for each design, using the
    # alignment produced above
    for i in $(seq 1 $NUM_DESIGNS); do
        if [ ! -f $outdir/designs/design-${i}.naive-design.tsv ]; then
            echo "design_naively.py $outdir/input-alns/design-${i}.fasta $outdir/designs/design-${i}.naive-design.tsv --window-size $ARG_WINDOWSIZE -gl $ARG_GL -gm $ARG_GM --verbose &> $outdir/designs/design-${i}.naive-design.out" >> $commands_fn
        fi
    done
}


# Run for Zika virus
run_for_taxid "64320" "None" "NC_035889,NC_012532" "AY632535"

# Run for Lassa virus, S segment
run_for_taxid "11620" "S" "KM821998,GU481072,KM821773" "KM821998"

# Run for Lassa virus, L segment
run_for_taxid "11620" "L" "U73034" "U73034"

# Run for Ebola virus (Zaire)
#run_for_taxid "186538" "None" "NC_002549"

# Run for Nipah virus
run_for_taxid "121791" "None" "NC_002728" "AF212302"

# Run for HIV-1
run_for_taxid "11676" "None" "NC_001802" "AF033819"

# Run for HCV
run_for_taxid "11103" "None" "NC_004102,NC_030791,NC_009827,NC_009826,NC_009825,NC_038882,NC_009824,NC_009823" "AF011751"

# Run for IAV segment 2
#run_for_taxid "11320" "2" "NC_026435,NC_002021,NC_007375,NC_026423,NC_007372" "GQ323558"

# Run for human coronavirus 229E
#run_for_taxid "11137" "None" "NC_028752"

# Run for Rhinovirus A, B, C
run_for_taxid "147711" "None" "NC_038311,NC_001617,NC_038311" "FJ445111"
run_for_taxid "147712" "None" "NC_038312" "DQ473485"
run_for_taxid "463676" "None" "NC_038878" "EF077279"

# Run for Enterovirus A, B, C, D
run_for_taxid "138948" "None" "NC_038306,NC_001612,NC_038306" "AY421760"
run_for_taxid "138949" "None" "NC_001472,NC_038307,NC_038307,NC_001472" "M88483"
run_for_taxid "138950" "None" "NC_002058" "V01149"
run_for_taxid "138951" "None" "NC_038308,NC_001430" "D00820"


# Load environment, and variables, for ADAPT
source ~/misc-repos/adapt-designs/scripts/run-adapt/custom-env/load_custom_env.sh

# Run parallel on the commands
echo "Running all commands.."
parallel --jobs $NJOBS --no-notice < $commands_fn
rm $commands_fn
