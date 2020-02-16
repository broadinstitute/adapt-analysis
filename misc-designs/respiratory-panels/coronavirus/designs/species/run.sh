#!/bin/bash

# Design panel.

# Make FASTA for all SARS (incl. SARS-CoV-2)
mkdir -p inputs
SARS_FASTA="inputs/all-sars.fasta"
source ~/anaconda3/etc/profile.d/conda.sh
conda activate adapt
python /ebs/dgd-analysis/dgd-analysis/utils/download_fasta_for_accessions.py ../../inputs/sars-minus-sarscov2.acc.txt > $SARS_FASTA
conda deactivate
zcat ../../inputs/sarscov2.20200214.fasta.gz >> $SARS_FASTA

FASTA_TO_USE="inputs/all-sars.manual-fasta.tsv"
echo -e "694009\tNone\tinputs/all-sars.fasta" > $FASTA_TO_USE

source ../../../run_common.sh

# Set input
IN_TSV="taxa.tsv"

ARG_IDM="4"
ARG_IDFRAC="0.01"

# Run the design
/usr/bin/time -f "mem=%K RSS=%M elapsed=%e cpu.sys=%S .user=%U" design.py complete-targets auto-from-file $IN_TSV $OUT_TSV_DIR -gl $ARG_GL -gm $ARG_GM -gp $ARG_GP -pl $ARG_PL -pm $ARG_PM -pp $ARG_PP --id-m $ARG_IDM --id-frac $ARG_IDFRAC --id-method shard --require-flanking3 H --max-primers-at-site $ARG_MAXPRIMERSATSITE --primer-gc-content-bounds 0.4 0.6 --max-target-length $ARG_MAXTARGETLENGTH --cost-fn-weights $ARG_COSTFNWEIGHTS --best-n-targets $ARG_BESTNTARGETS --predict-activity-model-path $PREDICTIVE_MODEL --mafft-path $MAFFT_PATH --prep-memoize-dir $PREP_MEMOIZE_DIR --cluster-threshold $CLUSTER_THRESHOLD --use-fasta $FASTA_TO_USE --write-input-seqs --write-input-aln --ncbi-api-key $NCBI_API_KEY --verbose &> out/design.out

# Re-run only SARS, with more strict parameters
ONLY_DESIGN_FOR="inputs/only-sars.tsv"
echo -e "694009\tNone" > $ONLY_DESIGN_FOR
ARG_GM="0"
ARG_PM="0"
ARG_GP="1.0"
ARG_IDM="5"
ARG_IDFRAC="0"
/usr/bin/time -f "mem=%K RSS=%M elapsed=%e cpu.sys=%S .user=%U" design.py complete-targets auto-from-file $IN_TSV $OUT_TSV_DIR -gl $ARG_GL -gm $ARG_GM -gp $ARG_GP -pl $ARG_PL -pm $ARG_PM -pp $ARG_PP --id-m $ARG_IDM --id-frac $ARG_IDFRAC --id-method shard --require-flanking3 H --max-primers-at-site $ARG_MAXPRIMERSATSITE --primer-gc-content-bounds 0.4 0.6 --max-target-length $ARG_MAXTARGETLENGTH --cost-fn-weights $ARG_COSTFNWEIGHTS --best-n-targets $ARG_BESTNTARGETS --predict-activity-model-path $PREDICTIVE_MODEL --mafft-path $MAFFT_PATH --prep-memoize-dir $PREP_MEMOIZE_DIR --cluster-threshold $CLUSTER_THRESHOLD --use-fasta $FASTA_TO_USE --only-design-for $ONLY_DESIGN_FOR --write-input-seqs --write-input-aln --ncbi-api-key $NCBI_API_KEY --verbose &> out/design.sars-strict.out

# gzip the stdout/stderr
gzip out/design.out
gzip out/design.sars-strict.out
