#!/bin/bash

# Args:
#   1: subcommand; 'index' or 'align'

# Allow activating conda environments
source ~/anaconda3/etc/profile.d/conda.sh


if [[ $1 == "index" ]]; then
    # Index the reference sequences to align against

    # Activate environment with bowtie2
    conda activate misc-bioinformatic-tools

    # Build an index of human transcripts (GENCODE v38) and common bacterial
    # pathogens
    # Also include, as positive controls, SARS-CoV-1 and SARS-CoV-2; these
    # also ensure there can be multiple alignments
    bowtie2-build data/gencode.v38.transcripts.fa.gz,data/bacterial-pathogens.fasta.gz,data/sars-cov-1-and-2.refseq.fasta.gz data/bowtie-index/background
elif [[ $1 == "align" ]]; then
    # Align experimentally tested guides to the background index

    # Activate environment with bowtie2
    conda activate misc-bioinformatic-tools

    bowtie2 -a -x data/bowtie-index/background -U data/experimentally-tested-guides.fasta -f -S out/aln-guides-to-background.sam --end-to-end -N 1 -L 7 -i S,1,1 --ma 0 --mp 1,1 --rdg 100,1 --rfg 100,1 --score-min L,-4,0 --threads 64
    # Below are reasons for various argument values:
    #   -a to report all alignments; helpful if index contains diversity of a pathogen, but very slow
    #   -N 1 -L 7 to pigeonhole 4 mismatches across the guide (at most 1 mismatch per seed), so we always find hit with 4 or fewer mismatches
    #   --end-to-end so whole guide aligns
    #   -i S,1,1 to stride seeds by 1+1*sqrt(28)=6 nt, i.e., so they basically do not overlap (which should not be needed with `-N 1 -L 7`)
    #   --ma 0 to not use a bonus score for matches (should be 0 anyway because --end-to-end is set)
    #   --mp 1,1 for penalty of 1 with mismatch
    #   --rdg 100,1 --rfg 100,1 for high penalty (at least 100) for gaps, i.e., disallow gaps in the alignment to match the behavior of ADAPT
    #   --score-min L,-4,0 to report alignments with score >= -4 (greatest possible alignment score, for a perfect guide-target match, is 0)
    # 
    # Runtime is very slow because of the value of -L; here are some runtime values
    # when tested against only human transcripts:
    #   -L 7: 12 hours
    #   -L 8: 58 mins
    #   -L 9: 5 mins
    #   -L 10: 40 sec

    # Compress the sam; a bam would be the normal way to go, but .sam.gz
    # requires less fancy tools to work with
    gzip -f out/aln-guides-to-background.sam
else
    echo "Unknown subcommand '$1'"
    exit 1
fi
