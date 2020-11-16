#!/bin/bash

# Summarize runtimes across windows using output of run.sh

for tax in sars-related-cov rhinovirus-a lassa-S; do
    SUMMARY=summary/summary.${tax}.tsv

    # Write header
    echo -e "run\tmemoize\ttime_from_start\twindow_num\twindow_start\twindow_end" > $SUMMARY

    for i in $(seq 1 3); do
        for memoize in yes no; do
            # Use the stdout file to determine runtimes
            outf="out/designs.${tax}_${i}_${memoize}.out.gz"

            # Find the timestamp of the start of the search (first window)
            # in milliseconds
            startdate=$(zcat $outf | grep "Found window" | head -n 1 | awk -F' - ' '{print $1}')
            starttime=$(date -d "$startdate" +%s%3N)

            windownum=1
            while read -r line; do
                # Calculate the time difference from the start in milliseconds
                datestr=$(echo "$line" | awk -F' - ' '{print $1}')
                timestamp=$(date -d "$datestr" +%s%3N)
                timediff=$(($timestamp-$starttime))

                # Find window
                window=$(echo "$line" | grep -oP 'Found window \[([0-9]+, [0-9]+)\)' | sed -e 's/Found window \[//' | sed -e 's/)//')
                windowstart=$(echo "$window " | awk -F', ' '{print $1}')
                windowend=$(echo "$window " | awk -F', ' '{print $2}')

                echo -e "$i\t$memoize\t$timediff\t$windownum\t$windowstart\t$windowend" >> $SUMMARY

                # Increment counter
                windownum=$((windownum+1))

                # Only process the first 100,000 windows
                if [[ "$windownum" -gt "100000" ]]; then
                    break
                fi
            done < <(zcat $outf | grep "Found window")
        done
    done

    gzip $SUMMARY
done
