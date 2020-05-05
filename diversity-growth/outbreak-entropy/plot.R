# Plot change in entropy over time, across genomic windows.
#
# This plots the genome on the x-axis and time on the y-axis, going
# downward.
#
# By Hayden Metsky <hayden@mit.edu>


require(ggplot2)
require(viridis)
require(ggpubr)

args <- commandArgs(trailingOnly=TRUE)
IN.TSV <- args[1]
OUT.PDF <- args[2]


# Read TSV of results and replace '_' in column names with '.'
entropy <- read.table(gzfile(IN.TSV), header=TRUE, sep="\t")
names(entropy) <- gsub("_", ".", names(entropy))

# Make midpoint value for each window
entropy$window.mid <- entropy$window.start + (entropy$window.end - entropy$window.start)/2

# Cutoff far left and right of genome, which have high entropy mostly
# due to technical artifacts -- i.e., drop the first and last window
entropy <- entropy[entropy$window.start != min(entropy$window.start),]
entropy <- entropy[entropy$window.start != max(entropy$window.start),]

# Make a hard cutoff on entropy of 0.3 to avoid issues from outliers
#entropy$entropy[entropy$entropy > 0.3] <- 0.3

# Make a heatmap
p <- ggplot(entropy, aes(y=start.date, x=window.mid, fill=entropy)) +
        geom_tile() +
        theme_pubr() +
        scale_y_discrete(limits = rev(levels(entropy$start.date))) +
        scale_fill_gradient(low="white", high="#3E1451") +
        xlab("Genome") +
        ylab("Date of sample collection")

ggsave(OUT.PDF, p, device="pdf", width=16, height=8, useDingbats=FALSE)
