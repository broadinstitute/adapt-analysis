#!/bin/bash

# Compute coverage of designs against sequences from each year.
#
# Author: Hayden Metsky <hayden@mit.edu>


# Allow activating conda environments
source ~/anaconda3/etc/profile.d/conda.sh

conda activate adapt

# Set variables for design
NJOBS=8
ARG_GM="1"
ARG_PM="3"

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
    commands_fn="/tmp/commands-covg-${taxid}_${segment}"
    echo -n "" > $commands_fn

    # Find the first year, N, with at least 10 sequences, and
    # compute coverage against sequences from all years <= N-1
    # (this ensures we can show something for sequences from years in START_YEAR..N-1;
    # otherwise there may be a lot of information we don't show if N
    # is a relatively recent year)
    for year_against in $(seq $START_YEAR $END_YEAR); do
        # Find the number of accessions for this year
        nseqs=$(cat $outdir/accessions.all-years.tsv | awk -v year="$year_against" '$4==year {print $3}' | wc -l)
        if [ "$nseqs" -ge "10" ]; then
            first_year_with_enough_seqs=$year_against
            break
        fi
    done
    initial_against="upto"`expr $first_year_with_enough_seqs - 1`

    # Compute coverage against the year $initial_against, as well as
    # each of the years in the range [$START_YEAR, $END_YEAR]
    years_against=("$initial_against" $(seq $START_YEAR $END_YEAR))

    # Produce commands that compute coverage, against accessions collected
    # in each year, for each design
    for year_against in "${years_against[@]}"; do
        # Compute coverage of each design against $year_against

        if [[ $year_against == upto* ]]; then
            # Pull out accessions up to the given year
            year_upto=$(echo "$year_against" | sed 's/upto//')
            cat $outdir/accessions.all-years.tsv | awk -v year="$year_upto" '$4<=year {print $3}' > $outdir/designs/accessions.in-${year_against}.txt
        else
            # Pull out accessions for this year
            cat $outdir/accessions.all-years.tsv | awk -v year="$year_against" '$4==year {print $3}' > $outdir/designs/accessions.in-${year_against}.txt
        fi
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
run_for_taxid "64320" "None" "NC_035889,NC_012532"

# Run for Lassa virus, S segment
run_for_taxid "11620" "S" "KM821998,GU481072,KM821773"

# Run for Ebola virus (Zaire)
run_for_taxid "186538" "None" "NC_002549"

# Run for Nipah virus
run_for_taxid "121791" "None" "NC_002728"

# Run for HIV-1
# Skip - this is very memory intensive
#run_for_taxid "11676" "None" "NC_001802"

# Run for HCV
run_for_taxid "11103" "None" "NC_004102,NC_030791,NC_009827,NC_009826,NC_009825,NC_038882,NC_009824,NC_009823"
