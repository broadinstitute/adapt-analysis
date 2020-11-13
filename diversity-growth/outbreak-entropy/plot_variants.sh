#!/bin/bash

# Allow activating conda environments
source ~/anaconda3/etc/profile.d/conda.sh

# Activate data-analysis
conda activate data-analysis

# There are a couple of ways to plot the data:
#   * Using cumulative genomes or not. Cumulative means that we count variants
#     using all genomes collected up to each date X, whereas not being cumulative means
#     only counting variants in genomes collected in the time window [X minus 7 days, X].
#   * Subsampling a number of genomes for each date or using all available genomes. This
#     is important to control for increased number of sequenced genomes -- without
#     accounting for this, the number of variants will of course be higher simply because
#     there are many more recent than older genomes.
# I opt to use cumulative genomes because (a) for some time intervals there are relatively few
# genomes, and variant counts may be noisy in these intervals or affected by large sequencing
# dumps, and (b) we may still care about detecting variation present earlier in the outbreak,
# so we care about variants in these genomes. I also opt to *not* subsample (i.e., use all
# genomes) and to control for increase numbers of genomes by only calling variants that are
# present in >=1% of genomes (done in plot_variants.R).

Rscript plot_variants.R sars-cov-2/variants.no-subsample.cumulative.tsv.gz sars-cov-2/variants.pdf &> sars-cov-2/variants.plot.out
