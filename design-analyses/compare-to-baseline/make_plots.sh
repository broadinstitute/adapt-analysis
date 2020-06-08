#!/bin/bash

# Make plots.


# Allow activating conda environments
source ~/anaconda3/etc/profile.d/conda.sh

conda activate data-analysis

function plot_for_taxid() {
    taxid="$1"
    segment="$2"

    Rscript plot.R tax-${taxid}_${segment} plots
}

plot_for_taxid "64320" "None"
plot_for_taxid "11620" "S"
plot_for_taxid "11620" "L"
plot_for_taxid "121791" "None"
plot_for_taxid "11676" "None"
plot_for_taxid "11103" "None"
plot_for_taxid "11320" "2"
plot_for_taxid "147711" "None"
plot_for_taxid "147712" "None"
plot_for_taxid "463676" "None"
plot_for_taxid "138948" "None"
plot_for_taxid "138949" "None"
plot_for_taxid "138950" "None"
plot_for_taxid "138951" "None"
