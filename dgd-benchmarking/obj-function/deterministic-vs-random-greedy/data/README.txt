To avoid having to download from NCBI and align, and to ensure I can subsample
but have the same genomes for each run, I did the following: for each taxon, I
ran ADAPT --write-input-aln to write the alignment. Those alignments are here.
This ensures the sequences are subsampled, curated, etc. and aligned the same
way.
  - For sars-related-cov, this used `auto-from-args 694009 None NC_045512,NC_004718`,
    and also `--subsample-seqs 250` to limit the input
  - For rhinovirus-a, this used `auto-from-args 147711 None NC_038311,NC_001617`
  - For lassa-S, this used `11620 S KM821998,GU481072,KM821773`

Example command:
  design.py sliding-window auto-from-args 11620 S KM821998,GU481072,KM821773 /dev/null -gl 28 --mafft-path ~/viral-ngs/viral-ngs-etc/conda-env/bin/mafft --write-input-aln lassa-S.fasta --verbose
