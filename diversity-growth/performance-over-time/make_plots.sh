#!/bin/bash

# Produce plots from designs.

# Allow activating conda environments
source ~/anaconda3/etc/profile.d/conda.sh

# Activate data-analysis
conda activate data-analysis

for subtype in H1Nany H3Nany HanyN1 HanyN2; do
    Rscript plot.R iav-${subtype}/designs.tsv plots/iav-${subtype}.pdf
done
