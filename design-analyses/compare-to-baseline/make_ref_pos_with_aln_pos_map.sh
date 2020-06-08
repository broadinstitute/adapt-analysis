#!/bin/bash

# Make map of reference position to alignment positions.

# Load environment, and variables, for ADAPT
source ~/misc-repos/adapt-designs/scripts/run-adapt/custom-env/load_custom_env.sh


function make_for_taxid() {
    taxid="$1"
    segment="$2"
    refacc="$3"

    echo "Making map for taxid $taxid (segment: $segment)" > /dev/tty

    # Check directory and alignment exists
    outdir="tax-${taxid}_${segment}"
    if [ ! -f $outdir/input-alns/all-accessions.fasta ]; then
        echo "FATAL: Cannot find alignment"
        exit 1
    fi

    python ../scripts/map_ref_pos_with_aln_pos.py $outdir/input-alns/all-accessions.fasta $refacc $outdir/input-alns/ref-pos-map.tsv
}

make_for_taxid "64320" "None" "AY632535"
make_for_taxid "11620" "S" "KM821998"
make_for_taxid "11620" "L" "U73034"
make_for_taxid "121791" "None" "AF212302"
make_for_taxid "11676" "None" "AF033819"
make_for_taxid "11103" "None" "AF011751"
make_for_taxid "11320" "2" "GQ323558"
make_for_taxid "147711" "None" "FJ445111"
make_for_taxid "147712" "None" "DQ473485"
make_for_taxid "463676" "None" "EF077279"
make_for_taxid "138948" "None" "AY421760"
make_for_taxid "138949" "None" "M88483"
make_for_taxid "138950" "None" "V01149"
make_for_taxid "138951" "None" "D00820"
