#!/bin/bash

for aln in $(find alignments -name '*.fasta'); do
    aln=$(echo "$aln" | sed 's/alignments\///' | sed 's/\.fasta//')

    # Calculate the fraction of sequences with missing data at each position
    python calc_frac_missing_over_genome.py alignments/${aln}.fasta plots/${aln}.tsv

    # Plot a distribution
    Rscript plot_frac_missing_over_genome.R plots/${aln}.tsv plots/${aln}
done
