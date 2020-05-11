#!/bin/bash
#
# Make table showing performance of k-mers selected from sequences
# in different years, against sequences from future years.

# Allow activating conda environments
source ~/anaconda3/etc/profile.d/conda.sh

# Activate adapt
conda activate adapt

for subtype in H1Nany H3Nany HanyN1 HanyN2; do
    for design_year in $(seq 2007 2019); do
        python -u pick_kmers_and_compute_coverage.py --in-fasta iav-${subtype}/in/iav-${subtype}.2005-present.aligned.fasta.gz --out-tsv iav-${subtype}/by-year/designs.${design_year}.tsv --mismatches 1 --design-year $design_year --years-included-in-design 3 &> iav-${subtype}/by-year/designs.${design_year}.out &
    done
done
