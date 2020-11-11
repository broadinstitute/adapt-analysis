#!/bin/bash

# Run with or without differential identification, for H or N subtypes.
#
# Args:
#  1: subtype to design for: 'H' or 'N'
#  2: whether or not to perform differential identification when designing
#     guides: 'id' or 'no-id'
subtype=$1
id=$2

function join { local IFS="$1"; shift; echo "$*"; }

H_subtypes=("H1" "H2" "H3" "H4" "H5" "H6" "H7" "H8" "H9" "H10" "H11" "H12" "H13" "H14" "H15" "H16")
N_subtypes=("N1" "N2" "N3" "N4" "N5" "N6" "N7" "N8" "N9")

# Unzip all input sequences
gzip -d data/iav-seqs/*.fasta.gz

# Make input
if [[ $subtype == "H" ]]; then
    for h in "${H_subtypes[@]}"; do
        in_fastas+=("data/iav-seqs/2k8_${h}Nstar_aligned_amplicon.fasta")
    done
    # The H amplicon length (only one amplicon is given as input) is 141 nt; use
    # a window of 111 nt so that there's a region in the middle where guides are
    # designed, away from the ends (which contain the primers)
    window_length="111"
elif [[ $subtype == "N" ]]; then
    for n in "${N_subtypes[@]}"; do
        in_fastas+=("data/iav-seqs/2k8_Hstar${n}_aligned_amplicon.fasta")
    done
    # The N amplicon length (only one amplicon is given as input) is 67 nt; use
    # a window of 37 nt so that there's a region in the middle where guides are
    # designed, away from the ends (which contain the primers)
    window_length="37"
else
    echo "Unknown subtype $subtype"
    exit 1
fi
input=$(join " " ${in_fastas[@]})

# Make output
if [[ $subtype == "H" ]]; then
    if [[ "$id" == "id" ]]; then
        for h in "${H_subtypes[@]}"; do
            out_tsvs+=("guide-designs/id/2k8_${h}Nstar.id.tsv")
        done
    elif [[ "$id" == "no-id" ]]; then
        for h in "${H_subtypes[@]}"; do
            out_tsvs+=("guide-designs/no-id/2k8_${h}Nstar.no-id.tsv")
        done
    else
        echo "Unknown id $id"
        exit 1
    fi
elif [[ $subtype == "N" ]]; then
    if [[ "$id" == "id" ]]; then
        for n in "${N_subtypes[@]}"; do
            out_tsvs+=("guide-designs/id/2k8_Hstar${n}.id.tsv")
        done
    elif [[ "$id" == "no-id" ]]; then
        for n in "${N_subtypes[@]}"; do
            out_tsvs+=("guide-designs/no-id/2k8_Hstar${n}.no-id.tsv")
        done
    else
        echo "Unknown id $id"
        exit 1
    fi
else
    echo "Unknown subtype $subtype"
fi
output=$(join " " ${out_tsvs[@]})
out_file="out/${subtype}.${id}.out"

# Construct a command to run, appending stderr and stdout to out_file
if [[ "$id" == "id" ]]; then
    cmd="python -u ~/diagnostic-guide-design/bin/design_guides.py $input -o $output -l 28 -w $window_length -m 1 -p 0.9 --id --id-m 3 --id-frac 0.05 --verbose >>$out_file 2>&1"
elif [[ "$id" == "no-id" ]]; then
    cmd="python -u ~/diagnostic-guide-design/bin/design_guides.py $input -o $output -l 28 -w $window_length -m 1 -p 0.9 --verbose >>$out_file 2>&1"
else
    echo "Unknown id $id"
    exit 1
fi

# Write command to out_file
echo -e "RUNNING: \"$cmd\"\n" > $out_file

# Add time and memory report to command
cmd="/usr/bin/time -f 'real %e\nuser %U\nsys %S\nmrss %M' $cmd"

# Run the command
eval $cmd
