#!/bin/bash

# adapt conda environment should be activated before running
# Args:
#   1: name of taxon (e.g., 'lasv-L')

accessions=$(cat data/$1/aln.fasta | grep '>' | sed 's/>//' | tr '\n' ' ')

python ../../utils/fetch_metadata.py -a $accessions -o out/$1/per-seq-metadata.tsv

