#!/bin/bash

# Args:
#   1: version (e.g., 'v3') to process
if [ -z "$1" ]; then
    echo "Version is not given"
    exit 1
fi
version=$1

mkdir -p tmp/${version}
../../../scripts/calculate_subtype_coverage.sh N data/${version}-guides/2k8${version}_N_crRNA_targets.fasta data/iav-seqs tmp/${version} out/${version}-guide-covg
