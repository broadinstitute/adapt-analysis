#!/bin/bash

# Make plots.


# Allow activating conda environments
source ~/anaconda3/etc/profile.d/conda.sh

conda activate data-analysis

function plot_for_taxid() {
    taxid="$1"
    segment="$2"
    tax_name="$3"

    Rscript plot.R tax-${taxid}_${segment} "$tax_name" plots
}

plot_for_taxid "64320" "None" "Zika virus"
plot_for_taxid "11620" "S" "Lassa virus, segment S"
plot_for_taxid "11620" "L" "Lassa virus, segment L"
plot_for_taxid "121791" "None" "Nipah virus"
plot_for_taxid "11676" "None" "HIV-1"
plot_for_taxid "11103" "None" "Hepacivirus C"
#plot_for_taxid "11320" "2" "Influenza A virus, segment 2"
plot_for_taxid "147711" "None" "Rhinovirus A"
plot_for_taxid "147712" "None" "Rhinovirus B"
plot_for_taxid "463676" "None" "Rhinovirus C"
plot_for_taxid "138948" "None" "Enterovirus A"
plot_for_taxid "138949" "None" "Enterovirus B"
plot_for_taxid "138950" "None" "Enterovirus C"
plot_for_taxid "138951" "None" "Enterovirus D"
