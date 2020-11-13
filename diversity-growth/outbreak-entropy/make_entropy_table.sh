#!/bin/bash
#
# Make table of entropy values for different time intervals
# over genomic windows.

# Allow activating conda environments
source ~/anaconda3/etc/profile.d/conda.sh

# Activate adapt
conda activate adapt

python compute_entropy_over_time.py --in-fasta sars-cov-2/in/gisaid_msa_2020-11-12/msa_1112.fasta.gz --out-tsv sars-cov-2/entropy.tsv --ref-acc 'EPI_ISL_402119' --start-date 2020-02-01 --end-date 2020-11-01 --interval-days 7 --window-length 500 --window-stride 25 -k 500 --sample-size-per-date-interval -1

gzip sars-cov-2/entropy.tsv
