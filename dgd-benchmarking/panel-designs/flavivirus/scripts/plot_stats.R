#!/usr/bin/env Rscript

# Make plots from summary stats of a design.
#
# This uses the summary stat output in out/
#
# By Hayden Metsky <hayden@mit.edu>

library(ggplot2)
library(gridExtra)
library(viridis)
library(ggforce)

summary.fn <- "out/summary.tsv"
summary <- data.frame(read.table(summary.fn, header=TRUE, sep="\t"))

# Replace '_' in column names with '.'
names(summary) <- gsub("_", ".", names(summary))

# Pull out taxa with a comparison to the non-speific design
summary.ns <- summary[!is.na(summary$nonspecific.cost), ]

# Compute differences between this design and the non-specific design
summary.ns$num.primers.total <- summary.ns$num.primers5 + summary.ns$num.primers3
summary.ns$nonspecific.num.primers.total <- summary.ns$nonspecific.num.primers5 + summary.ns$nonspecific.num.primers3
summary.ns$cost.diff <- summary.ns$cost - summary.ns$nonspecific.cost
summary.ns$target.len.diff <- summary.ns$target.len - summary.ns$nonspecific.target.len
summary.ns$num.primers.total.diff <- summary.ns$num.primers.total - summary.ns$nonspecific.num.primers.total
summary.ns$num.guides.diff <- summary.ns$num.guides - summary.ns$nonspecific.num.guides

# Add dummy x-axis variable for below
summary.ns$x.const <- 1
summary.ns$x.const <- factor(summary.ns$x.const)

# Plot difference in cost
p <- ggplot(summary.ns, aes(x=x.const, y=cost.diff))
p <- p + geom_jitter(width=0.4, height=0, size=5)
p <- p + ylab("Difference in cost") + ylim(-0.5, 1.0)
p <- p + scale_x_discrete(breaks=NULL)  # remove x-axis; it's meaningless here
p <- p + theme_bw(base_size=18) # bw & larger font sizes
p <- p + theme(panel.grid.minor=element_blank(), # leave out minor grid
               aspect.ratio=2,  # vertically stretched
               axis.title.x=element_blank())    # no x-axis title
ggsave("plots/cost-diff.pdf", p, width=8, height=8, useDingbats=FALSE)

# Plot difference in target length
p <- ggplot(summary.ns, aes(x=x.const, y=target.len.diff))
p <- p + geom_jitter(width=0.4, height=0, size=5)
p <- p + ylab("Difference in target length (nt)") + ylim(-5, 30)
p <- p + scale_x_discrete(breaks=NULL)  # remove x-axis; it's meaningless here
p <- p + theme_bw(base_size=18) # bw & larger font sizes
p <- p + theme(panel.grid.minor=element_blank(), # leave out minor grid
               aspect.ratio=2,  # vertically stretched
               axis.title.x=element_blank())    # no x-axis title
ggsave("plots/target-len-diff.pdf", p, width=8, height=8, useDingbats=FALSE)

# Plot difference in total number of primers
p <- ggplot(summary.ns, aes(x=x.const, y=num.primers.total.diff))
p <- p + geom_jitter(width=0.4, height=0, size=5)
p <- p + ylab("Difference in total number of primers") + ylim(-1, 2)
p <- p + scale_x_discrete(breaks=NULL)  # remove x-axis; it's meaningless here
p <- p + theme_bw(base_size=18) # bw & larger font sizes
p <- p + theme(panel.grid.minor=element_blank(), # leave out minor grid
               aspect.ratio=2,  # vertically stretched
               axis.title.x=element_blank())    # no x-axis title
ggsave("plots/total-num-primers-diff.pdf", p, width=8, height=8, useDingbats=FALSE)

# Plot difference in number of guides
p <- ggplot(summary.ns, aes(x=x.const, y=num.guides.diff))
p <- p + geom_jitter(width=0.4, height=0, size=5)
p <- p + ylab("Difference in number of guides") + ylim(-1, 2)
p <- p + scale_x_discrete(breaks=NULL)  # remove x-axis; it's meaningless here
p <- p + theme_bw(base_size=18) # bw & larger font sizes
p <- p + theme(panel.grid.minor=element_blank(), # leave out minor grid
               aspect.ratio=2,  # vertically stretched
               axis.title.x=element_blank())    # no x-axis title
ggsave("plots/num-guides-diff.pdf", p, width=8, height=8, useDingbats=FALSE)
