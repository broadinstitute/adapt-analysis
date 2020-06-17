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

summary.fn <- "summary/summary.tsv"
summary <- data.frame(read.table(summary.fn, header=TRUE, sep="\t"))

# Replace '_' in column names with '.'
names(summary) <- gsub("_", ".", names(summary))

summary$name <- factor(summary$name)

# Pull out taxa with a comparison to the non-speific design
summary.ns <- summary[!is.na(summary$nonspecific.cost), ]

# Compute differences between this design and the non-specific design
summary.ns$num.primers.total <- summary.ns$num.primers5 + summary.ns$num.primers3
summary.ns$nonspecific.num.primers.total <- summary.ns$nonspecific.num.primers5 + summary.ns$nonspecific.num.primers3
summary.ns$cost.diff <- summary.ns$cost - summary.ns$nonspecific.cost
summary.ns$target.len.diff <- summary.ns$target.len - summary.ns$nonspecific.target.len
summary.ns$num.primers.total.diff <- summary.ns$num.primers.total - summary.ns$nonspecific.num.primers.total
summary.ns$num.guides.diff <- summary.ns$num.guides - summary.ns$nonspecific.num.guides

# Determine non-specific and specific colors
nonspecific.pt.color <- viridis(5)[1]
specific.pt.color <- viridis(5)[4]

# The dengue cost is so much larger than the others and throws off the
# cost plot; bring down values > 3 by 2
summary.ns$cost <- ifelse(summary.ns$cost<=3, summary.ns$cost, ifelse(summary.ns$cost<5,NA,summary.ns$cost-2))
summary.ns$nonspecific.cost <- ifelse(summary.ns$nonspecific.cost<=3, summary.ns$nonspecific.cost, ifelse(summary.ns$nonspecific.cost<5,NA,summary.ns$nonspecific.cost-2))

# Plot non-specific and specific costs
p <- ggplot(summary.ns, aes(x=name))
p <- p + geom_point(aes(y=nonspecific.cost), color=nonspecific.pt.color, size=4)
p <- p + geom_point(aes(y=cost), color=specific.pt.color, size=4)
p <- p + ylab("Cost") + xlab("Species")
p <- p + scale_y_continuous(breaks=c(2.5,3,3.5,4), labels=c(2.5,3,5.5,6), limits=c(2.5,4))  # accommodate high dengue cost value by inserting break
p <- p + theme_bw(base_size=18) # bw & larger font sizes
p <- p + theme(panel.grid.minor=element_blank(), # leave out minor grid
               aspect.ratio=1,  # square
               axis.text.x=element_text(angle=45, vjust=1, hjust=1)) # 45 degree x-axis label
ggsave("plots/costs.pdf", p, width=8, height=8, useDingbats=FALSE)

# Plot non-specific and specific target lengths
p <- ggplot(summary.ns, aes(x=name))
p <- p + geom_point(aes(y=nonspecific.target.len), color=nonspecific.pt.color, size=4)
p <- p + geom_point(aes(y=target.len), color=specific.pt.color, size=4)
p <- p + scale_y_continuous(breaks=c(80,90,100,110,120), limits=c(80,120))
p <- p + ylab("Target length") + xlab("Species")
p <- p + theme_bw(base_size=18) # bw & larger font sizes
p <- p + theme(panel.grid.minor=element_blank(), # leave out minor grid
               aspect.ratio=1,  # square
               axis.text.x=element_text(angle=45, vjust=1, hjust=1)) # 45 degree x-axis label
ggsave("plots/target-len.pdf", p, width=8, height=8, useDingbats=FALSE)

# Plot non-specific and specific total number of primers
p <- ggplot(summary.ns, aes(x=name))
p <- p + geom_point(aes(y=nonspecific.num.primers.total), color=nonspecific.pt.color, size=4)
p <- p + geom_point(aes(y=num.primers.total), color=specific.pt.color, size=4)
p <- p + ylab("Total number of primers") + xlab("Species")
p <- p + theme_bw(base_size=18) # bw & larger font sizes
p <- p + theme(panel.grid.minor=element_blank(), # leave out minor grid
               aspect.ratio=1,  # square
               axis.text.x=element_text(angle=45, vjust=1, hjust=1)) # 45 degree x-axis label
ggsave("plots/total-num-primers.pdf", p, width=8, height=8, useDingbats=FALSE)

# Plot non-specific and specific number of guides
p <- ggplot(summary.ns, aes(x=name))
p <- p + geom_point(aes(y=nonspecific.num.guides), color=nonspecific.pt.color, size=4)
p <- p + geom_point(aes(y=num.guides), color=specific.pt.color, size=4)
p <- p + scale_y_continuous(breaks=c(1,2,3,4,5,6))
p <- p + ylab("Number of guides") + xlab("Species")
p <- p + theme_bw(base_size=18) # bw & larger font sizes
p <- p + theme(panel.grid.minor=element_blank(), # leave out minor grid
               aspect.ratio=1,  # square
               axis.text.x=element_text(angle=45, vjust=1, hjust=1)) # 45 degree x-axis label
ggsave("plots/num-guides.pdf", p, width=8, height=8, useDingbats=FALSE)
