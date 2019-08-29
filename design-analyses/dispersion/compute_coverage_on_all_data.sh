#!/bin/bash

# Compute coverage of designs against all sequences.
#
# The designs using resampled data are, of course, designed against only a subset
# of all the data. This computes the coverage that these designs achieve against
# all of the sequences.
#
# Author: Hayden Metsky <hayden@mit.edu>


# Allow activating conda environments
source ~/anaconda3/etc/profile.d/conda.sh

conda activate adapt

# Set variables for computing coverage
NJOBS=72
ARG_GM="1"
ARG_PM="3"


function run_for_taxid() {
    # Set information on the taxonomy, from arguments
    taxid="$1"
    segment="$2"

    echo "Computing coverage for $taxid (segment: $segment)" > /dev/tty

    outdir="tax-${taxid}_${segment}"

    # Make sure there exists accessions for this taxid
    if [ ! -f $outdir/accessions.tsv ]; then
        echo "FATAL: Nonexistent accessions.tsv file in $outdir"
        exit 1
    fi

    # Pull out only the accessions
    cat $outdir/accessions.tsv | awk '{print $3}' > $outdir/accessions.acc-only.txt.tmp

    # Make a directory in which to place the coverage data
    mkdir -p $outdir/designs/resampled/coverages

    # Write commands to a file
    commands_fn="/tmp/commands-dispersioncovg-${taxid}_${segment}"
    echo -n "" > $commands_fn

    # For every design, compute coverage against all the accessions
    for design_name in $(ls -1 $outdir/designs/resampled/ | grep '.tsv.0' | sed 's/\.tsv.0//'); do
        echo "analyze_coverage.py $outdir/designs/resampled/${design_name}.tsv.0 $outdir/accessions.acc-only.txt.tmp -gm $ARG_GM -pm $ARG_PM --use-accessions --write-frac-bound $outdir/designs/resampled/coverages/${design_name}.coverage-against-all.txt --verbose &> $outdir/designs/resampled/coverages/${design_name}.coverage-against-all.out" >> $commands_fn
    done

    # Run parallel
    parallel --jobs $NJOBS --no-notice < $commands_fn

    rm $commands_fn
    rm $outdir/accessions.acc-only.txt.tmp
}


# Run for Zika virus
run_for_taxid "64320" "None"

# Run for Lassa virus, S segment
run_for_taxid "11620" "S"

# Run for Ebola virus (Zaire)
run_for_taxid "186538" "None" "NC_002549"

# Run for Nipah virus
run_for_taxid "121791" "None" "NC_002728"

# Run for HIV-1
run_for_taxid "11676" "None" "NC_001802"

# Run for HCV
# See make_designs.sh for why this is skipped
#run_for_taxid "11103" "None" "NC_004102,NC_030791,NC_009827,NC_009826,NC_009825,NC_038882,NC_009824,NC_009823"
