I previously produced alignments via ~/tmp/ncov-gisaid/align.sh, which uses
MAFFT with the command:
```
mafft --preservecase --auto --keeplength --addfragments [input].fasta NC_045512.2.fasta > [output].fasta
```

For the Nov. 12 2020 download from GISAID, it seems GISAID capped the number
of sequences I could download to 10,000. So, instead, I downloaded the MSA
that they provided directly. That is in gisaid_msa_2020-11-12/. This has
180,318 sequences from which the scripts can parse a date (and thus include
in annalyses); 3,879 sequences did not have a parseable date.
