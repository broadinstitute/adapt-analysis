#!/bin/bash

# Compile a distribution of coverages by year for each taxonomy.
#
# Author: Hayden Metsky <hayden@mit.edu>


# Let the coverage of a design be the mean of the top TOPN targets
# in the design
TOPN=20

# Set variables for measuring change over time
START_YEAR=2005
END_YEAR=2019


for dir in $(ls -1 . | grep 'tax-'); do
    design_name=$(echo "$dir" | sed -s 's/tax-//')
    echo "Compiling coverages for $design_name"

    # Start a file with the distribution of coverage on the test set
    echo -e "coverage.against.year\tnum.seqs.against.year\tdesign.up.to.year\tnum.seqs.up.to.year\tcoverage" > tax-${design_name}/coverages.distribution.txt

    # Find all the years against which coverage was computed
    years_against=( $(ls -1 tax-${design_name}/designs/accessions.in-*.txt | xargs -n1 basename | awk -F'in-' '{print $2}' | sed 's/\.txt//' | sort -n) )

    for year_against in "${years_against[@]}"; do
        # Find coverage of each design against sequences collected in $year_against

        if [ ! -s tax-${design_name}/designs/accessions.in-${year_against}.txt ]; then
            # There are no accessions collected in $year_against
            # Skip this year
            continue
        fi

        for year_design in $(seq $START_YEAR $END_YEAR); do
            # Find coverage of each design produced using sequences up to
            # and including $year_design, against sequences collected in
            # $year_against

            # Determine number of sequences up to year $year_design and at year $year_against
            num_seqs_year_against=$(cat tax-${design_name}/designs/accessions.in-${year_against}.txt | wc -l)
            num_seqs_year_design=$(cat tax-${design_name}/designs/accessions.up-to-${year_design}.tsv | wc -l)

            # Compute a coverage value for each resampling
            for design_file in $(ls -1 tax-${design_name}/designs/designs_up-to-${year_design}/coverages/coverage-in-${year_against}/design-*.coverage.txt); do
                # Take mean of the coverage values of the $TOPN targets
                topn_mean=$(cat $design_file | tail -n +2 | awk -F '\t' -v topn="$TOPN" '$1<=topn {s+=$5; n+=1} END {print s/n}')
                echo -e "${year_against}\t${num_seqs_year_against}\t${year_design}\t${num_seqs_year_design}\t${topn_mean}" >> tax-${design_name}/coverages.distribution.txt
            done
        done
    done
done
