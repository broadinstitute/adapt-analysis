#!/bin/bash

# Print summary statistics for the Lassa segment S example
# Each statistic is a mean across the 5 samplings of sequences

DESIGNS="tax-11620_S/designs"

echo -n "ADAPT, 1 probe, mean over genome: "
for i in $(seq 1 5); do paste $DESIGNS/design-$i.real-design.maximize-activity.hgc-1.tsv $DESIGNS/design-$i.naive-design.mode-upto-1.tsv | tail -n +2 | awk -F'\t' '{print $5}' | awk '{s+=$1}END{print s/NR}';done | awk '{s+=$1}END{print s/NR}'

echo -n "ADAPT, 1 probe, median over genome: "
for i in $(seq 1 5); do paste $DESIGNS/design-$i.real-design.maximize-activity.hgc-1.tsv $DESIGNS/design-$i.naive-design.mode-upto-1.tsv | tail -n +2 | awk -F'\t' '{print $5}' | sort -g | awk ' { a[i++]=$1; } END { print a[int(i/2)]; }';done | awk '{s+=$1}END{print s/NR}'

echo -n "Mode, mean over genome: "
for i in $(seq 1 5); do paste $DESIGNS/design-$i.real-design.maximize-activity.hgc-1.tsv $DESIGNS/design-$i.naive-design.mode-upto-1.tsv | tail -n +2 | awk -F'\t' '{print $17}' | awk '{s+=$1}END{print s/NR}';done | awk '{s+=$1}END{print s/NR}'

echo -n "Mode, median over genome: "
for i in $(seq 1 5); do paste $DESIGNS/design-$i.real-design.maximize-activity.hgc-1.tsv $DESIGNS/design-$i.naive-design.mode-upto-1.tsv | tail -n +2 | awk -F'\t' '{print $17}' | sort -g | awk ' { a[i++]=$1; } END { print a[int(i/2)]; }';done | awk '{s+=$1}END{print s/NR}'

echo -n "Consensus, mean over genome: "
for i in $(seq 1 5); do paste $DESIGNS/design-$i.real-design.maximize-activity.hgc-1.tsv $DESIGNS/design-$i.naive-design.consensus.tsv | tail -n +2 | awk -F'\t' '{print $16}' | awk '{s+=$1}END{print s/NR}';done | awk '{s+=$1}END{print s/NR}'

echo -n "Consensus, median over genome: "
for i in $(seq 1 5); do paste $DESIGNS/design-$i.real-design.maximize-activity.hgc-1.tsv $DESIGNS/design-$i.naive-design.consensus.tsv | tail -n +2 | awk -F'\t' '{print $16}' | sort -g | awk ' { a[i++]=$1; } END { print a[int(i/2)]; }';done | awk '{s+=$1}END{print s/NR}'

echo -n "ADAPT 1 probe - mode, mean over genome: "
for i in $(seq 1 5); do paste $DESIGNS/design-$i.real-design.maximize-activity.hgc-1.tsv $DESIGNS/design-$i.naive-design.mode-upto-1.tsv | tail -n +2 | awk -F'\t' '{print $5-$17}' | awk '{s+=$1}END{print s/NR}';done | awk '{s+=$1}END{print s/NR}'

echo -n "ADAPT 1 probe - mode, median over genome: "
for i in $(seq 1 5); do paste $DESIGNS/design-$i.real-design.maximize-activity.hgc-1.tsv $DESIGNS/design-$i.naive-design.mode-upto-1.tsv | tail -n +2 | awk -F'\t' '{print $5-$17}' | sort -g | awk ' { a[i++]=$1; } END { print a[int(i/2)]; }';done | awk '{s+=$1}END{print s/NR}'

echo -n "ADAPT 1 probe - consensus, mean over genome: "
for i in $(seq 1 5); do paste $DESIGNS/design-$i.real-design.maximize-activity.hgc-1.tsv $DESIGNS/design-$i.naive-design.consensus.tsv | tail -n +2 | awk -F'\t' '{print $5-$16}' | awk '{s+=$1}END{print s/NR}';done | awk '{s+=$1}END{print s/NR}'

echo -n "ADAPT 1 probe - consensus, median over genome: "
for i in $(seq 1 5); do paste $DESIGNS/design-$i.real-design.maximize-activity.hgc-1.tsv $DESIGNS/design-$i.naive-design.consensus.tsv | tail -n +2 | awk -F'\t' '{print $5-$16}' | sort -g | awk ' { a[i++]=$1; } END { print a[int(i/2)]; }';done | awk '{s+=$1}END{print s/NR}'
