#!/bin/bash
#
# Make table of number of variants at each genome position at different
# time points.

# Allow activating conda environments
source ~/anaconda3/etc/profile.d/conda.sh

# Activate adapt
conda activate adapt

python compute_variants_over_time.py --in-fasta sars-cov-2/in/gisaid_msa_2020-11-12/msa_1112.fasta.gz --out-tsv sars-cov-2/variants.no-subsample.no-cumulative.tsv --ref-acc 'EPI_ISL_402119' --start-date 2020-02-01 --end-date 2020-11-01 --interval-by-month --sample-size-per-date-interval -1
python compute_variants_over_time.py --in-fasta sars-cov-2/in/gisaid_msa_2020-11-12/msa_1112.fasta.gz --out-tsv sars-cov-2/variants.no-subsample.cumulative.tsv --ref-acc 'EPI_ISL_402119' --start-date 2020-02-01 --end-date 2020-11-01 --interval-by-month --sample-size-per-date-interval -1 --cumulative
python compute_variants_over_time.py --in-fasta sars-cov-2/in/gisaid_msa_2020-11-12/msa_1112.fasta.gz --out-tsv sars-cov-2/variants.subsample-300.no-cumulative.tsv --ref-acc 'EPI_ISL_402119' --start-date 2020-02-01 --end-date 2020-11-01 --interval-by-month --sample-size-per-date-interval 300
python compute_variants_over_time.py --in-fasta sars-cov-2/in/gisaid_msa_2020-11-12/msa_1112.fasta.gz --out-tsv sars-cov-2/variants.subsample-300.cumulative.tsv --ref-acc 'EPI_ISL_402119' --start-date 2020-02-01 --end-date 2020-11-01 --interval-by-month --sample-size-per-date-interval 300 --cumulative

gzip sars-cov-2/variants.no-subsample.no-cumulative.tsv
gzip sars-cov-2/variants.no-subsample.cumulative.tsv
gzip sars-cov-2/variants.subsample-300.no-cumulative.tsv
gzip sars-cov-2/variants.subsample-300.cumulative.tsv
