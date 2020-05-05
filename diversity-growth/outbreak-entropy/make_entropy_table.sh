#!/bin/bash
#
# Make table of entropy values for different time intervals
# over genomic windows.

# Allow activating conda environments
source ~/anaconda3/etc/profile.d/conda.sh

# Activate adapt
conda activate adapt

python compute_entropy_over_time.py --in-fasta sars-cov-2/in/gisaid_hcov-19_2020_04_28_06.aligned.fasta.gz --out-tsv sars-cov-2/entropy.tsv --ref-acc 'EPI_ISL_402119' --start-date 2020-03-01 --end-date 2020-04-07 --interval-days 7 --window-length 500 --window-stride 500 --sample-size-per-date-interval 500
