#!/usr/bin/env Rscript

# Plot a distribution over positions in an alignment of the
# fraction of sequences that have missing data.
#
# Args:
#  1: TSV file; col 1 gives position in an alignment and
#     col 2 gives fraction of sequences with missing data
#  2: prefix of output PDF file; outputs .hist.pdf (histogram)
#     and .density.pdf (KDE)
#
# By Hayden Metsky <hayden@mit.edu>

library(ggplot2)

# Parse arguments
args <- commandArgs(trailingOnly=TRUE)
in.tsv <- args[1]
out.pdf.prefix <- args[2]

# Read input data
data <- read.table(in.tsv, header=TRUE)

# Plot a histogram
p <- ggplot(data, aes(frac.missing))
p <- p + geom_histogram()
p <- p + xlab("Fraction of sequences with missing data") +
	ylab("Number of positions")

# Leave out usual ggplot2 background and grid lines but keep border
# Use aspect.ratio=1 to make the plot square
simple_theme <- theme_bw() +
		 theme(panel.grid.major=element_blank(),
               panel.grid.minor=element_blank(),
               strip.background=element_blank(),
               panel.border=element_rect(colour="black"),
               aspect.ratio=1)
p <- p + simple_theme

p <- p + ggsave(paste0(out.pdf.prefix, ".hist.pdf"),
	width=8, height=8, useDingbats=FALSE)

# Plot a KDE
p <- ggplot(data, aes(frac.missing))
p <- p + geom_density(fill="gray", adjust=2)
p <- p + xlim(0, 1)
p <- p + xlab("Fraction of sequences with missing data") +
	ylab("Density")
p <- p + simple_theme
p <- p + ggsave(paste0(out.pdf.prefix, ".density.pdf"),
	width=8, height=8, useDingbats=FALSE)