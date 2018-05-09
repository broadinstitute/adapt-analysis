#!/bin/bash

# Args:
#  1: subtype (H or N)
#  2: input guides (fasta)
#  3: path to directory containing fasta files of sequences against which to calculate
#     coverage (e.g., if subtype is H, this directory should contain the files
#     2k8_[SUBTYPE]Nstar.fasta where [SUBTYPE] is H1, H2, etc.)
#  4: path to tmp directory
#  5: path to output directory

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

guides=$2
seqsdir=$3
tmpdir=$4
outdir=$5

# Make a fastq of the guide sequences
echo -n "" > $tmpdir/guides.fastq
for line in $(cat $guides); do
    if [[ $line == \>* ]]; then
        # Line is a sequence header
        header=$(echo "$line" | sed 's/>/@/')
        echo "$header" >> $tmpdir/guides.fastq
    else
        # Line is a sequence
        echo "$line" >> $tmpdir/guides.fastq
        echo "+" >> $tmpdir/guides.fastq

        # Print artificial quality score (~ times the
        # guide length, 28)
        echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~" >> $tmpdir/guides.fastq
    fi
done

# Put all guide names into a file
grep '>' $guides | sed 's/>//' > $tmpdir/guide-names.txt

# Find the coverage that the guides give across each subtype
# individually
for subtype in "${SUBTYPES[@]}"; do
    # Use bwa to align each guide against the sequences for this
    # subtype
    # The parameters require a total of 21 bases to match for an alignment,
    # with no tolerance for gaps
    # The options '-k 8 -r 0.0001 -c 1000000 -s 1000000 -y 1000000' are needed for
    # bwa to find all alignments when the reference genomes (in $seqsdir/*.fasta) are highly
    # similar; e.g., without these, it will toss MEMs that occur too many times in the
    # genomes (which might be necessary as seeds), will filter chains that would be
    # necessary to keep, etc. The option '-y' (for third reseeding) slows this down
    # considerably, but its value needs to be high to find/keep seeds that occur often
    # in the (highly similar) target genomes.
    if [[ $1 == "H" ]]; then
        reffasta="$seqsdir/2k8_${subtype}Nstar.fasta"
    elif [[ $1 == "N" ]]; then
        reffasta="$seqsdir/2k8_Hstar${subtype}.fasta"
    else
        echo "Unknown subtype $1"
        exit 1
    fi
    bwa mem -a -k 8 -r 0.0001 -c 1000000 -s 1000000 -y 1000000 -M -A 1 -B 0 -O 1000000 -E 1000000 -L 1000000 -T 21 $reffasta $tmpdir/guides.fastq > $tmpdir/guides-to-${subtype}.all.sam

    # Only keep mapped guides
    samtools view -F 4 $tmpdir/guides-to-${subtype}.all.sam > $tmpdir/guides-to-${subtype}.sam
    rm $tmpdir/guides-to-${subtype}.all.sam

    num_target_seqs=$(grep '>' $reffasta | wc -l)

    # Summarize results for different thresholds on the number of
    # total matching bases in the alignment (here: alignment score, or AS)
    # Since '-T 21' was passed to bwa above, using as=21 should use all output alignments
    for as in 21 23 25 27; do
        python $(dirname "$0")/summarize_guide_mapping.py $tmpdir/guides-to-${subtype}.sam $tmpdir/guide-names.txt $num_target_seqs $tmpdir/guides-to-${subtype}.summary.as${as}.per-guide.tsv.tmp $tmpdir/guides-to-${subtype}.summary.as${as}.per-guide-set.tsv.tmp --aln-score-filter $as

        # Add a column to each output TSV (in the middle) giving the subtype
        awk -v subtype="$subtype" '{print $1"\t"subtype"\t"$2}' $tmpdir/guides-to-${subtype}.summary.as${as}.per-guide.tsv.tmp > $tmpdir/guides-to-${subtype}.summary.as${as}.per-guide.tsv
        rm $tmpdir/guides-to-${subtype}.summary.as${as}.per-guide.tsv.tmp
        awk -v subtype="$subtype" '{print $1"\t"subtype"\t"$2}' $tmpdir/guides-to-${subtype}.summary.as${as}.per-guide-set.tsv.tmp > $tmpdir/guides-to-${subtype}.summary.as${as}.per-guide-set.tsv
        rm $tmpdir/guides-to-${subtype}.summary.as${as}.per-guide-set.tsv.tmp
    done
done

# Combine the output TSVs across subtypes into one large table
for as in 21 23 25 27; do
    # Combine fractions per-guide and per-guide-set
    for k in per-guide per-guide-set; do
        # Make file header
        outfn="$outdir/covg.as${as}.${k}.tsv"
        echo -n "name" > $outfn
        for subtype in "${SUBTYPES[@]}"; do
            echo -ne "\t$subtype" >> $outfn
        done
        echo "" >> $outfn

        cat $tmpdir/guides-to-*.summary.as${as}.${k}.tsv > $tmpdir/guides-to-ALL.summary.as${as}.${k}.tsv
        for name in $(cat $tmpdir/guides-to-ALL.summary.as${as}.${k}.tsv | awk '{print $1}' | sort | uniq); do
            # Make a separate file with just this guide or guide set (the one with name)
            awk -v name="$name" '$1==name {print $0}' $tmpdir/guides-to-ALL.summary.as${as}.${k}.tsv > $tmpdir/guides-to-ALL.summary.as${as}.${k}.name-${name}.tsv

            # Write a row for this guide or guide set (with each column being a subtype)
            echo -n "$name" >> $outfn
            for subtype in "${SUBTYPES[@]}"; do
                frac=$(awk -v subtype="$subtype" '$2==subtype {print $3}' $tmpdir/guides-to-ALL.summary.as${as}.${k}.name-${name}.tsv)
                echo -ne "\t$frac" >> $outfn
            done
            echo "" >> $outfn
            
            rm $tmpdir/guides-to-ALL.summary.as${as}.${k}.name-${name}.tsv
        done
    done
done
