#!/bin/bash

# Make designs to use for comparing performance against baseline (naive) designs.
#
# To make a more clear comparison to naive designs (which design one guide per
# window), this only uses dgd with the sliding window approach.
#
# Author: Hayden Metsky <hayden@mit.edu>


# Allow activating conda environments
source ~/anaconda3/etc/profile.d/conda.sh

# Set variables for measuring uncertainty
NUM_DESIGNS=10

# Set variables for design
NJOBS=2
PREP_MEMOIZE_DIR="/ebs/dgd-analysis/prep-memoize-dir"
MAFFT_PATH="/home/hayden/viral-ngs/viral-ngs-etc/conda-env/bin/mafft"
CLUSTER_THRESHOLD=1.0   # Use high value to obtain a single cluster
ARG_GL="28"
ARG_GM="1"
ARG_GP="0.99"
ARG_WINDOWSIZE="200"


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
    mkdir -p $outdir/input-alns
    mkdir -p $outdir/designs

    # Fetch accessions and create table of them
    conda activate data-analysis
    if [ ! -f $outdir/accessions.tsv ]; then
        python ../scripts/find_year_for_accessions.py $taxid $segment | awk -v taxid="$taxid" -v segment="$segment" '{print taxid"\t"segment"\t"$1}' | sort | uniq > $outdir/accessions.tsv
    fi

    # Activate environment for design
    conda activate dgd

    # Write commands to a file
    commands_fn="/tmp/commands-designs-${taxid}_${segment}"
    echo -n "" > $commands_fn

    # Randomly sample a number of sequences equal to the number
    # of accessions (bootstrapping over the input)
    sample_size=$(cat $outdir/accessions.tsv | wc -l)

    # Produce a design.py command for each design (each with randomly sampled input)
    # Be sure to save the exact input (alignment) it used, to use this as input
    # later with the naive approach
    for i in $(seq 1 $NUM_DESIGNS); do
        echo "design.py sliding-window auto-from-args $taxid $segment $refaccs $outdir/designs/design-${i}.real-design.tsv --window-size $ARG_WINDOWSIZE -gl $ARG_GL -gm $ARG_GM -gp $ARG_GP --mafft-path $MAFFT_PATH --prep-memoize-dir $PREP_MEMOIZE_DIR --sample-seqs $sample_size --cluster-threshold $CLUSTER_THRESHOLD --use-accessions $outdir/accessions.tsv --write-input-aln $outdir/input-alns/design-${i}.fasta --verbose &> $outdir/designs/design-${i}.real-design.out" >> $commands_fn
    done

    # Run parallel on the design.py commands
    parallel --jobs $NJOBS --no-notice < $commands_fn

    echo -n "" > $commands_fn

    # Produce a design_naively.py command for each design, using the
    # same input that was used with design.py (the input is
    # design-${i}.fasta.0, with the `.0` because there will be only
    # one cluster)
    for i in $(seq 1 $NUM_DESIGNS); do
        echo "design_naively.py $outdir/input-alns/design-${i}.fasta.0 $outdir/designs/design-${i}.naive-design.tsv --window-size $ARG_WINDOWSIZE -gl $ARG_GL -gm $ARG_GM --verbose &> $outdir/designs/design-${i}.naive-design.out" >> $commands_fn
    done

    # Run parallel on the design_naively.py commands
    parallel --jobs $NJOBS --no-notice < $commands_fn

    rm $commands_fn
}


# Run for Zika virus
#run_for_taxid "64320" "None" "NC_035889,NC_012532"

# Run for Lassa virus, S segment
#run_for_taxid "11620" "S" "KM821998,GU481072,KM821773"

# Run for Hepatitis C virus (Hepacivirus C)
run_for_taxid "11103" "None" "NC_004102,NC_030791,NC_009827,NC_009826,NC_009825,NC_038882,NC_009824,NC_009823"
