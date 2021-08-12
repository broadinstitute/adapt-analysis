#!/bin/bash

# Make plots.


# Allow activating conda environments
source ~/anaconda3/etc/profile.d/conda.sh

conda activate data-analysis

Rscript plot.R standard plots/cross-validation.standard.all-taxonomies.pdf
Rscript plot.R relaxed plots/cross-validation.relaxed.all-taxonomies.pdf
rm -f Rplots.pdf
