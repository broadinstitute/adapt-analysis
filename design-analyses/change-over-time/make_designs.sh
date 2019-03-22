#!/bin/bash

# Make designs to use for seeing how coverage changes over time.
#
# Author: Hayden Metsky <hayden@mit.edu>


# Allow activating conda environments
source ~/anaconda3/etc/profile.d/conda.sh

# Make tmp directory for memoizing alignments and stats
mkdir -p /tmp/prep-memoize-dir/

# Set variables for measuring uncertainty
NUM_DESIGNS=100

# Set variables for design
NJOBS=16
MAFFT_PATH="/home/hayden/viral-ngs/viral-ngs-etc/conda-env/bin/mafft"
CLUSTER_THRESHOLD=1.0   # Use high value to obtain a single cluster
ARG_GL="28"
ARG_GM="1"
ARG_GP="0.99"
ARG_PL="30"
ARG_PM="3"
ARG_PP="0.95"
ARG_MAXPRIMERSATSITE="5"
ARG_MAXTARGETLENGTH="1000"
ARG_COSTFNWEIGHTS="0.6667 0.2222 0.1111"
ARG_BESTNTARGETS="20"

# Set variables for measuring change over time
MIN_YEAR=1900
START_YEAR=2005
END_YEAR=2019


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
    if [ ! -f $outdir/accessions.all-years.tsv ]; then
        python ../scripts/find_year_for_accessions.py $taxid $segment | awk -v taxid="$taxid" -v segment="$segment" '{print taxid"\t"segment"\t"$1"\t"$2}' | sort | uniq > $outdir/accessions.all-years.tsv
    fi

    # Activate environment for design
    conda activate dgd

    # Write commands to a file
    commands_fn="/tmp/commands-${taxid}_${segment}"
    echo -n "" > $commands_fn

    # Produce commands that design using all accessions collected up to each year
    for year in $(seq $START_YEAR $END_YEAR); do
        # Make an output directory for designs for this year
        mkdir -p $outdir/designs/designs_up-to-${year}

        # Pull out all accessions collected up to and including this year
        cat $outdir/accessions.all-years.tsv | awk -v min="$MIN_YEAR" -v max="$year" '$4 >= min && $4 <= max {print $1"\t"$2"\t"$3}' > $outdir/designs/accessions.up-to-${year}.tsv

        # Randomly sample a number of sequences equal to the number
        # of accessions (bootstrapping over the input)
        sample_size=$(cat $outdir/designs/accessions.up-to-${year}.tsv | wc -l)

        # Produce a command for each design (each with randomly sampled input)
        for i in $(seq 1 $NUM_DESIGNS); do
            echo "design.py complete-targets auto-from-args $taxid $segment $refaccs $outdir/designs/designs_up-to-${year}/design-${i}.tsv -gl $ARG_GL -gm $ARG_GM -gp $ARG_GP -pl $ARG_PL -pm $ARG_PM -pp $ARG_PP --max-primers-at-site $ARG_MAXPRIMERSATSITE --max-target-length $ARG_MAXTARGETLENGTH --cost-fn-weights $ARG_COSTFNWEIGHTS --best-n-targets $ARG_BESTNTARGETS --mafft-path $MAFFT_PATH --prep-memoize-dir /tmp/prep-memoize-dir --sample-seqs $sample_size --cluster-threshold $CLUSTER_THRESHOLD --use-accessions $outdir/designs/accessions.up-to-${year}.tsv --verbose &> $outdir/designs/designs_up-to-${year}/design-${i}.out" >> $commands_fn
        done
    done

    # Run parallel
    parallel --jobs $NJOBS --no-notice < $commands_fn

    rm $commands_fn
}


# Run for Zika virus
run_for_taxid "64320" "None" "NC_035889,NC_012532"
