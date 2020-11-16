#!/bin/bash

# Produce designs.
# ADAPT must be loaded (conda env) before running.

ADAPT_PATH="$HOME/adapt"


# Make tmp file with commands
commandsfile=$(mktemp)
echo "Saving commands to tmp file $commandsfile"

echo "Making designs"

# Since there is randomness, do this all 3 times; note, however, that
# inputs are the same
# Try 3 taxa: sars-related-cov should work with 1-2 guides, and
# rhinovirus-a has a lot of diversity and may require many guides,
# as well as lassa-S
echo -n "" > $commandsfile
for i in $(seq 1 3); do
    for memoize in yes no; do
        for tax in sars-related-cov rhinovirus-a lassa-S; do
            outf="designs.${tax}_${i}_${memoize}"
            if [ ! -f "designs/$outf.tsv" ]; then
                if [ "$memoize" == "yes" ]; then
                    memoize_arg=""
                elif [ "$memoize" == "no" ]; then
                    memoize_arg="--do-not-memoize-guide-computations"
                fi
                cmd="/usr/bin/time -f 'mem=%K RSS=%M elapsed=%e cpu.sys=%S .user=%U' design.py complete-targets fasta data/${tax}.fasta.gz -o designs/$outf.tsv -pm 3 -pp 0.9 --primer-gc-content-bounds 0.3 0.7 --max-primers-at-site 10 -gl 28 --max-target-len 250 --best-n-targets 10 --predict-activity-model-path $ADAPT_PATH/models/classify/model-51373185 $ADAPT_PATH/models/regress/model-f8b6fd5d --obj maximize-activity --soft-guide-constraint 1 --hard-guide-constraint 5 --penalty-strength 0.25 --maximization-algorithm random-greedy --specific-against-taxa specificity-lists/${tax}.tsv --id-m 4 --id-frac 0.01 --id-method shard $memoize_arg --verbose &> out/$outf.out"
                echo "$cmd" >> $commandsfile
            fi
        done
    done
done

# Run commands in parallel
parallel --jobs 10 --no-notice < $commandsfile

rm $commandsfile

echo "Done making designs"

# Make a summary file that has the elapsed runtime for each design
#summaryfile="designs.summary.tsv"
#echo -e "run\ttaxonomy\tmemoize\telapsed_runtime" > $summaryfile
#for i in $(seq 1 3); do
#    for memoize in yes no; do
#        for tax in sars-related-cov rhinovirus-a lassa-S; do
#            outf="designs.${tax}_${i}_${memoize}"
#            runtime=$(tail -n 1 out/$outf.out | awk '{print $3}' | sed 's/elapsed=//')
#            echo "$i\t$tax\t$memoize\t$runtime" >> $summaryfile
#        done
#    done
#done
#gzip -f $summaryfile

gzip out/*.out
