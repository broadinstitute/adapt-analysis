#!/bin/bash

python ../download_seqs.py make-species-list taxid10239.nbr.20180905.modified.txt -o species-list.tsv --human-host-lineages-to-add human-host-lineages-to-add.txt --influenza-seqs human-influenza-ABC-seqs.tsv
