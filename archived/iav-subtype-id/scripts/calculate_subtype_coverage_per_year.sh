#!/bin/bash

# Args:
#  1: subtype (H or N)
#  2: path to directory containing fasta files of sequences against which to calculate
#     coverage (e.g., if subtype is H, this directory should contain the files
#     2k8_[SUBTYPE]Nstar.fasta where [SUBTYPE] is H1, H2, etc.)
#  3: path to tsv file mapping accession of sequence (col 1) to year of
#     sequence (col 2)
#  4: path to tmp directory
#  5: path to output directory
#
# This uses the output in tmp/ generated by calculate_subtype_coverage.sh

H_SUBTYPES=( "H1" "H2" "H3" "H4" "H5" "H6" "H7" "H8" "H9" "H10" "H11" "H12" "H13" "H14" "H15" "H16" )
N_SUBTYPES=( "N1" "N2" "N3" "N4" "N5" "N6" "N7" "N8" "N9" )
if [[ $1 == "H" ]]; then
    SUBTYPES=("${H_SUBTYPES[@]}")
elif [[ $1 == "N" ]]; then
    SUBTYPES=("${N_SUBTYPES[@]}")
else
    echo "Unknown subtype $1"
    exit 1
fi

seqsdir=$2
yearstsv=$3
tmpdir=$4
outdir=$5

for subtype in "${SUBTYPES[@]}"; do
    if [[ $1 == "H" ]]; then
        seqs="$seqsdir/2k8_${subtype}Nstar.fasta"
    elif [[ $1 == "N" ]]; then
        seqs="$seqsdir/2k8_Hstar${subtype}.fasta"
    else
        echo "Unknown subtype $1"
        exit 1
    fi

    # Pull out alignments with score AS >= 27
    cat $tmpdir/guides-to-${subtype}.sam | awk '{split($14, subfield, ":"); if(subfield[3]>=27) print $0}' > $tmpdir/guides-to-${subtype}.orig-score.as27.sam

    # Pull out accessions from this subtype
    grep '>' $seqs | awk '{print $1}' | sed 's/>//' | sort > $tmpdir/${subtype}.acc.txt

    # Make an array of all guide sets for which to write a table (plus 'all' to write one
    # that combines all guide sets)
    readarray -t guide_sets_to_process < <(cat $tmpdir/guide-names.txt | awk -F'-' '{print $1"-"$2}' | sort | uniq)
    guide_sets_to_process+=("ALL")

    for guide_set in "${guide_sets_to_process[@]}"; do
        if [ "$guide_set" == "ALL" ]; then
            inguides="$tmpdir/guides-to-${subtype}.orig-score.as27.sam"
            outtsv="$outdir/covg-by-year.orig-score.as27.${subtype}.tsv"
        else
            # Filter the SAM for only guides from this guide set
            inguides="$tmpdir/guides-to-${subtype}.orig-score.as27.${guide_set}.sam"
            pattern="^${guide_set}-"
            awk -v pattern="$pattern" '$1 ~ pattern' $tmpdir/guides-to-${subtype}.orig-score.as27.sam > $inguides

            outtsv="$outdir/covg-by-year.orig-score.as27.${subtype}.${guide_set}.tsv"
        fi
            

        # Go through each unique year
        echo -e "year\tnum.seqs\tnum.covered\tfrac.covered" > $outtsv
        while read year; do
            # Find all accessions from year
            cat $yearstsv | awk -v y="$year" '$2==y {print $1}' | sort > $tmpdir/${year}.acc.txt

            # Find only those accessions/years from this subtype
            join -t $'\t' -1 1 -2 1 <(sort $tmpdir/${subtype}.acc.txt) <(sort $tmpdir/${year}.acc.txt) > $tmpdir/${subtype}.${year}.acc.txt

            numseqs=$(wc -l $tmpdir/${subtype}.${year}.acc.txt | awk '{print $1}')
            numcovered=$(comm -12 $tmpdir/${subtype}.${year}.acc.txt <(cat $inguides | awk '{print $3}' | sort) | wc -l)
            if [ "$numseqs" -eq "0" ]; then
                # There are no sequences for this year in this subtype; none
                # are covered
                frac="0"
            else
                frac=$(echo "scale=4; $numcovered/$numseqs" | bc)
            fi

            echo -e "$year\t$numseqs\t$numcovered\t$frac" >> $outtsv

            rm $tmpdir/${subtype}.${year}.acc.txt
            rm $tmpdir/${year}.acc.txt
        done < <(cat $yearstsv | awk '{print $2}' | sed '/^$/d' | sort -n | uniq)

        if [ "$guide_set" != "ALL" ]; then
            rm $inguides
        fi
    done

    rm $tmpdir/guides-to-${subtype}.orig-score.as27.sam
    rm $tmpdir/${subtype}.acc.txt
done
