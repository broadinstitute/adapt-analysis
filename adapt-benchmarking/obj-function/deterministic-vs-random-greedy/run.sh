#!/bin/bash

# Produce designs.
# ADAPT must be loaded (conda env) before running.

ADAPT_PATH="$HOME/adapt"


# Make tmp file with commands
commandsfile=$(mktemp)
echo "Saving commands to tmp file $commandsfile"

echo "Making designs"

# Try different values for soft and hard guide constraints, and
# for the penalty strength
# Since there is randomness, do this all 5 times
# Try 3 taxa: sars-related-cov should work with 1-2 guides, and
# rhinovirus-a has a lot of diversity and may require many guides,
# as well as lassa-S
echo -n "" > $commandsfile
for i in $(seq 1 5); do
    for hgc in $(seq 1 4); do
        for ((sgc=1; sgc<=hgc; sgc++)); do
            for ps in 0.1 0.5 1.0; do
                for algo in greedy random-greedy; do
                    for tax in sars-related-cov rhinovirus-a lassa-S; do
                        outf="designs.${tax}_${algo}_${hgc}_${sgc}_${ps}_${i}.tsv"
                        if [ ! -f "designs/$outf" ]; then
                            cmd="design.py complete-targets fasta data/${tax}.fasta.gz -o designs/$outf -pm 3 -pp 0.9 --primer-gc-content-bounds 0.3 0.7 --max-primers-at-site 10 -gl 28 --max-target-len 250 --best-n-targets 10 --predict-activity-model-path $ADAPT_PATH/models/classify/model-51373185 $ADAPT_PATH/models/regress/model-f8b6fd5d --obj maximize-activity --soft-guide-constraint $sgc --hard-guide-constraint $hgc --penalty-strength $ps --maximization-algorithm $algo --verbose &> out/$outf.out"
                            echo "$cmd" >> $commandsfile
                        fi
                    done
                done
            done
        done
    done
done

# Run commands in parallel
parallel --jobs 24 --no-notice < $commandsfile

rm $commandsfile

echo "Done making designs"

# Make a summary file that has the top 5 design options for
# each choice of parameter values
summaryfile="designs.summary.tsv"
echo -e "run\thard_guide_constraint\tsoft_guide_constraint\tpenalty_strength\talgorithm\ttaxonomy\ttarget_objective_value\tnum_guides\texpected_activity" > $summaryfile
for i in $(seq 1 5); do
    for hgc in $(seq 1 4); do
        for ((sgc=1; sgc<=hgc; sgc++)); do
            for ps in 0.1 0.5 1.0; do
                for algo in greedy random-greedy; do
                    for tax in sars-related-cov rhinovirus-a lassa-S; do
                        designfile="designs/designs.${tax}_${algo}_${hgc}_${sgc}_${ps}_${i}.tsv"
                        # Pull out the columns we want for the first 5 rows (excluding header)
                        cat $designfile | tail -n +2 | head -n 5 | awk -F'\t' -v i="$i" -v hgc="$hgc" -v sgc="$sgc" -v ps="$ps" -v algo="$algo" -v tax="$tax" '{print i"\t"hgc"\t"sgc"\t"ps"\t"algo"\t"tax"\t"$1"\t"$13"\t"$15}' >> $summaryfile
                    done
                done
            done
        done
    done
done
gzip -f $summaryfile
