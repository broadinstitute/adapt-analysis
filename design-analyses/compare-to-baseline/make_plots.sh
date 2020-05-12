#!/bin/bash

# Make plots.


# Allow activating conda environments
source ~/anaconda3/etc/profile.d/conda.sh

conda activate data-analysis

function plot_for_taxid() {
    taxid="$1"
    segment="$2"

    Rscript plot.R tax-${taxid}_${segment} plots/tax-${taxid}_${segment}.pdf
}

plot_for_taxid "64320" "None"
plot_for_taxid "11620" "S"
