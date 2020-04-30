#!/bin/bash

python scripts/select_taxa_from_accession_list.py 2020-04/taxid10239.20200429.with-additions.nbr.gz 2020-04/taxa.unmodified.tsv --force-include-taxa 2020-04/taxa-to-force-include.tsv  --force-exclude-taxa 2020-04/taxa-to-exclude.tsv
