#!/bin/bash

# Design panel.

source ../../../run_common.sh

# Set input
IN_TSV="taxa.tsv"

ARG_IDM="4"
ARG_IDFRAC="0.01"

# Make accessions input file
ACCESSIONS_TO_USE="accessions-to-use.tsv"
cat ../../inputs/IAV.Hany-Nany.segment2.acc.txt | grep . | awk '{print "11320\t2\t"$1}' > $ACCESSIONS_TO_USE
cat ../../inputs/IBV.segment1.acc.txt | grep . | awk '{print "11520\t1\t"$1}' >> $ACCESSIONS_TO_USE

# Run the design
/usr/bin/time -f "mem=%K RSS=%M elapsed=%e cpu.sys=%S .user=%U" design.py complete-targets auto-from-file $IN_TSV $OUT_TSV_DIR -gl $ARG_GL -gm $ARG_GM -gp $ARG_GP -pl $ARG_PL -pm $ARG_PM -pp $ARG_PP --id-m $ARG_IDM --id-frac $ARG_IDFRAC --id-method shard --require-flanking3 H --max-primers-at-site $ARG_MAXPRIMERSATSITE --primer-gc-content-bounds 0.4 0.6 --max-target-length $ARG_MAXTARGETLENGTH --cost-fn-weights $ARG_COSTFNWEIGHTS --best-n-targets $ARG_BESTNTARGETS --predict-activity-model-path $PREDICTIVE_MODEL --mafft-path $MAFFT_PATH --prep-memoize-dir $PREP_MEMOIZE_DIR --cluster-threshold $CLUSTER_THRESHOLD --use-accessions $ACCESSIONS_TO_USE --write-input-seqs --write-input-aln --ncbi-api-key $NCBI_API_KEY --verbose &> out/design.out

# gzip the stdout/stderr
gzip out/design.out
