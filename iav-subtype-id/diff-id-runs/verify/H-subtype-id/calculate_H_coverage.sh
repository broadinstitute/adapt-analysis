#!/bin/bash

mkdir -p tmp/v3
../../../scripts/calculate_subtype_coverage.sh H data/v3-guides/2k8v3_H_crRNA_targets.fasta data/iav-seqs tmp/v3 out/v3-guide-covg
