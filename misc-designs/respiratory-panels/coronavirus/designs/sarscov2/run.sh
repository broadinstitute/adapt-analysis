#!/bin/bash

# Design panel.

source ../../../run_common.sh

# Note that the taxid 999999 is made-up, referring to (SARS except (SARS-CoV-2)).

# Be specific against all coronaviruses, but only design for sub-SARS lineages
ONLY_DESIGN_FOR="inputs/only-sars.tsv"

# Use accessions for the non-(SARS-CoV-2) strains
ACCESSIONS_TO_USE="inputs/sars-minus-sarscov2.tsv"
cat ../../inputs/sars-minus-sarscov2.acc.txt | grep . | awk '{print "999999\tNone\t"$1}' > $ACCESSIONS_TO_USE

# Use FASTA for SARS-CoV-2
FASTA_TO_USE="inputs/sarscov2.manual-fasta.tsv"
echo -e "2697049\tNone\t../../inputs/sarscov2.20200214.fasta.gz" > $FASTA_TO_USE

# Set input
IN_TSV="taxa.tsv"

# Set more strict parameters
ARG_GM="0"
ARG_PM="0"
ARG_GP="1.0"
ARG_IDM="5"
ARG_IDFRAC="0"

/usr/bin/time -f "mem=%K RSS=%M elapsed=%e cpu.sys=%S .user=%U" design.py complete-targets auto-from-file $IN_TSV $OUT_TSV_DIR -gl $ARG_GL -gm $ARG_GM -gp $ARG_GP -pl $ARG_PL -pm $ARG_PM -pp $ARG_PP --id-m $ARG_IDM --id-frac $ARG_IDFRAC --id-method shard --require-flanking3 H --max-primers-at-site $ARG_MAXPRIMERSATSITE --primer-gc-content-bounds 0.4 0.6 --max-target-length $ARG_MAXTARGETLENGTH --cost-fn-weights $ARG_COSTFNWEIGHTS --best-n-targets $ARG_BESTNTARGETS --predict-activity-model-path $PREDICTIVE_MODEL --mafft-path $MAFFT_PATH --prep-memoize-dir $PREP_MEMOIZE_DIR --cluster-threshold $CLUSTER_THRESHOLD --use-accessions $ACCESSIONS_TO_USE --use-fasta $FASTA_TO_USE --only-design-for $ONLY_DESIGN_FOR --write-input-seqs --write-input-aln --ncbi-api-key $NCBI_API_KEY --verbose &> out/design.out

# gzip the stdout/stderr
gzip out/design.out
