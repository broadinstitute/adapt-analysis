#!/bin/bash

# Make plots.


# Allow activating conda environments
source ~/anaconda3/etc/profile.d/conda.sh

conda activate data-analysis

Rscript plot.R plots/compare-to-baseline.all-taxonomies.pdf
