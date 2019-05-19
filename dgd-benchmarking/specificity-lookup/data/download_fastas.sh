#!/bin/bash

# Download a FASTA for each taxonomy.

# Allow activating conda environments
source ~/anaconda3/etc/profile.d/conda.sh

# download_fasta_for_taxonmy.py requires dgd
conda activate dgd

TAXONOMIES_FILE="taxonomies.tsv"

# Download for each taxonomy serially; they could be done
# in parallel but serially may be more kind to NCBI's servers
while read -r taxonomy; do
    taxid=$(echo "$taxonomy" | awk -F'\t' '{print $4}')
    segment=$(echo "$taxonomy" | awk -F'\t' '{print $5}')

    python ../../../utils/download_fasta_for_taxonomy.py $taxid $segment > fastas/${taxid}_${segment}.fasta
    gzip fastas/${taxid}_${segment}.fasta
done < <(tail -n +2 $TAXONOMIES_FILE)
