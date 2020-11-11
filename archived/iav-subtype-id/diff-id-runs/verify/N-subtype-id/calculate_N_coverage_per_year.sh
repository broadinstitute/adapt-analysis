#!/bin/bash

# Args:
#   1: version (e.g., 'v3') to process
if [ -z "$1" ]; then
    echo "Version is not given"
    exit 1
fi
version=$1

../../../scripts/calculate_subtype_coverage_per_year.sh N data/iav-seqs ../../data/iav-seqs-years/2k8_Hstar.acc-with-year.tsv tmp/${version} out/${version}-guide-covg
