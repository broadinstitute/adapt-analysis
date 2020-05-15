#!/bin/bash

#!/bin/bash

# Concatenate by-year TSV files into one TSV file for each subtype.


for subtype in H1Nany H3Nany HanyN1 HanyN2; do
    first_f=$(ls -1 iav-${subtype}/by-year/designs.*.tsv | head -n 1)
    # Concatenate TSV files, keeping only the header from the first
    (head -1 $first_f && tail -n +2 -q iav-${subtype}/by-year/designs.*.tsv) > iav-${subtype}/designs.tsv
done
