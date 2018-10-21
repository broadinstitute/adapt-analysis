#!/bin/bash

python ../download_seqs.py make-fasta-files taxid10239.nbr.20180905.modified.txt -o fasta-files --human-host-lineages-to-add human-host-lineages-to-add.txt --reverse-complement-accessions seqs-to-reverse-complement.txt --influenza-seqs human-influenza-ABC-seqs.tsv
