#!/bin/bash

# Make a plot for each taxonomy

for tax in sars-related-cov rhinovirus-a lassa-S; do
    Rscript plot.R summary/summary.${tax}.tsv.gz plots/plot.${tax}.pdf
done
