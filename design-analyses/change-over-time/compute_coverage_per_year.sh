#!/bin/bash

# Compute coverage of designs against sequences from each year.
#
# Author: Hayden Metsky <hayden@mit.edu>


# Allow activating conda environments
source ~/anaconda3/etc/profile.d/conda.sh

# Make tmp directory for memoizing alignments and stats
mkdir -p /tmp/prep-memoize-dir/

# Set variables for design
NJOBS=16
MAFFT_PATH="/home/hayden/viral-ngs/viral-ngs-etc/conda-env/bin/mafft"
CLUSTER_THRESHOLD=1.0   # Use high value to obtain a single cluster
ARG_GM="1"
ARG_PM="2"

# Set variables for measuring change over time
START_YEAR=2005
END_YEAR=2019


function run_for_taxid() {
    # Set information on the taxonomy, from arguments
    taxid="$1"
    segment="$2"

    echo "Computing coverage for $taxid (segment: $segment)" > /dev/tty

    outdir="tax-${taxid}_${segment}"

    # Make sure there exists accessions for this taxid
    if [ ! -f $outdir/accessions.all-years.tsv ]; then
        echo "FATAL: Nonexistent accessions.all-years.tsv file in $outdir"
        exit 1
    fi

    # Write commands to a file
    commands_fn="/tmp/commands-${taxid}_${segment}"
    echo -n "" > $commands_fn

    # Produce commands that compute coverage, against accessions collected
    # in each year, for each design
    for year_against in $(seq $START_YEAR $END_YEAR); do
        # Compute coverage of each design against $year_against

        # Pull out accessions for this year
        cat $outdir/accessions.all-years.tsv | awk -v year="$year_against" '$4==year {print $3}' > $outdir/designs/accessions.in-${year_against}.txt
        if [ ! -s $outdir/designs/accessions.in-${year_against}.txt ]; then
            # There are no accessions collected in $year_against
            continue
        fi

        # For every design, compute coverage against the accessions
        # collected in $year_against
        for year_design in $(seq $START_YEAR $END_YEAR); do
            mkdir -p $outdir/designs/designs_up-to-${year_design}/coverages/coverage-in-${year_against}

            for design_name in $(ls -1 $outdir/designs/designs_up-to-${year_design}/ | grep '.tsv.0' | sed 's/\.tsv.0//'); do
                echo "analyze_coverage.py $outdir/designs/designs_up-to-${year_design}/${design_name}.tsv.0 $outdir/designs/accessions.in-${year_against}.txt -gm $ARG_GM -pm $ARG_PM --use-accessions --write-frac-bound $outdir/designs/designs_up-to-${year_design}/coverages/coverage-in-${year_against}/${design_name}.coverage.txt --verbose &> $outdir/designs/designs_up-to-${year_design}/coverages/coverage-in-${year_against}/${design_name}.out" >> $commands_fn
            done
        done
    done

    # Run parallel
    parallel --jobs $NJOBS --no-notice < $commands_fn

    rm $commands_fn
}


# Run for Zika virus
run_for_taxid "64320" "None"

# Run for Lassa virus, S segment
run_for_taxid "11620" "S"
