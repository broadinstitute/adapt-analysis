#!/bin/bash

# Args:
#  1: filename listing a sequence on each line

# This assumes the year can be found in the pattern
#  '/[year] [4 digits]' where [year] is 4 digits and [4 digits] is any 4 digits

while read -r seqname; do
    # Print the sequence name in the first column
    echo -ne "$seqname\t"

    # Extract the year and print it in the second column (or 'unknown' if not found)
    # This will pull out [4 digits]_[4 digits] and only take the first [4 digits] as the
    # year
    # Use 'head -1' to only print the first match in a line (there should only be one match,
    # but in rare exceptions there may be more than one, in which case we may output an
    # incorrect year)
    (echo "$seqname" | grep -Po '/\d{4} \d{4}' || printf "unknown\n") | head -1 | sed 's/\///' | awk -F' ' '{print $1}'
done < "$1"
