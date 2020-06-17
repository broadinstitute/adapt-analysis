#!/usr/bin/env Rscript

# Make plots from summary stats of a design.
#
# This uses the summary stat output in data/
#
# By Hayden Metsky <hayden@mit.edu>

library(ggplot2)
library(gridExtra)
library(viridis)
library(ggforce)

summary.fn <- "data/design-results.tsv"
summary <- data.frame(read.table(summary.fn, header=TRUE, sep="\t"))

# Replace '_' in column names with '.'
names(summary) <- gsub("_", ".", names(summary))

# Add some stats
summary$frac.kept <- summary$num.curated.seqs / summary$num.input.seqs

# Add elapsed time and memory in new units
summary$elapsed.time.min <- summary$elapsed.time / 60.0
summary$rss.mb <- summary$rss / 1000.0

# Round mean number of guides, and convert to factor
summary$mean.cluster.num.guides.rounded <- round(summary$mean.cluster.num.guides)
summary[summary$mean.cluster.num.guides.rounded >= 9, ]$mean.cluster.num.guides.rounded <- "ge9"
summary$mean.cluster.num.guides.rounded <- factor(summary$mean.cluster.num.guides.rounded)

# Plot histogram of number of clusters per species
p <- ggplot(summary, aes(num.clusters))
p <- p + geom_histogram(binwidth=1)
p <- p + scale_x_continuous(breaks=seq(1, 15, by=2))
p <- p + xlab("Number of clusters") + ylab("Number of species")
p <- p + theme_bw(base_size=18) # bw & larger font sizes
p <- p + theme(panel.grid.minor=element_blank(), # leave out minor grid
               aspect.ratio=1)  # square
ggsave("plots/num-clusters.pdf", p, width=8, height=8, useDingbats=FALSE)

# Plot number of clusters vs. log(number of input sequences)
p <- ggplot(summary, aes(x=num.input.seqs, y=num.clusters))
p <- p + geom_point()
p <- p + scale_x_continuous(trans='log10')
p <- p + scale_y_continuous(breaks=seq(1, 5, by=1))
p <- p + xlab("Number of sequences") + ylab("Number of clusters")
p <- p + theme_bw(base_size=18) # bw & larger font sizes
p <- p + theme(panel.grid.minor=element_blank(), # leave out minor grid
               aspect.ratio=1)  # square
ggsave("plots/num-clusters-vs-num-seqs.pdf", p, width=8, height=8, useDingbats=FALSE)

# Plot histogram of the fraction of sequences kept after curation
p <- ggplot(summary, aes(frac.kept))
p <- p + geom_histogram(binwidth=0.05)
p <- p + xlab("Fraction of sequences kept") + ylab("Number of species")
p <- p + theme_bw(base_size=18) # bw & larger font sizes
p <- p + theme(panel.grid.minor=element_blank(), # leave out minor grid
               aspect.ratio=1)  # square
ggsave("plots/frac-sequences-kept.pdf", p, width=8, height=8, useDingbats=FALSE)

# Plot fraction of sequences kept vs. number of input sequences
p <- ggplot(summary, aes(x=num.input.seqs, y=frac.kept))
p <- p + geom_point()
p <- p + scale_x_continuous(trans='log10')
p <- p + xlab("Number of sequences") + ylab("Fraction of sequences kept")
p <- p + theme_bw(base_size=18) # bw & larger font sizes
p <- p + theme(panel.grid.minor=element_blank(), # leave out minor grid
               aspect.ratio=1)  # square
ggsave("plots/frac-sequences-kept-vs-num-seqs.pdf", p, width=8, height=8, useDingbats=FALSE)

# Plot mean genome region (amplicon) length vs. number of input seqs
# (here, means are taken over clusters)
p <- ggplot(summary, aes(x=num.input.seqs, y=mean.cluster.target.len))
p <- p + geom_point(aes(color=mean.cluster.num.guides))
p <- p + scale_x_continuous(trans='log10')
p <- p + xlab("Number of sequences") + ylab("Mean target length")
p <- p + scale_color_viridis(name="num.guides")  # viridis color scheme
p <- p + theme_bw(base_size=18) # bw & larger font sizes
p <- p + theme(panel.grid.minor=element_blank(), # leave out minor grid
               aspect.ratio=1)  # square
ggsave("plots/mean-cluster-target-len-vs-num-seqs.pdf", p, width=8, height=8, useDingbats=FALSE)

# Plot mean number of guides vs. number of input seqs
# (here, means are taken over clusters)
p <- ggplot(summary, aes(x=num.input.seqs, y=mean.cluster.num.guides))
p <- p + geom_point(aes(color=mean.cluster.target.len))
p <- p + scale_x_continuous(trans='log10')
p <- p + scale_y_continuous(breaks=seq(1, 20, by=5))
p <- p + xlab("Number of sequences") + ylab("Mean number of guides")
p <- p + scale_color_viridis(name="target.len")  # viridis color scheme
p <- p + theme_bw(base_size=18) # bw & larger font sizes
p <- p + theme(panel.grid.minor=element_blank(), # leave out minor grid
               aspect.ratio=1)  # square
ggsave("plots/mean-cluster-num-guides-vs-num-seqs.pdf", p, width=8, height=8, useDingbats=FALSE)

# Plot rounded mean number of guides vs. number of input seqs as sina plot
p <- ggplot(summary, aes(y=num.input.seqs, x=mean.cluster.num.guides.rounded))
p <- p + geom_sina(aes(group=mean.cluster.num.guides.rounded,
                       color=mean.cluster.target.len))
p <- p + scale_y_continuous(trans='log10')
p <- p + ylab("Number of sequences") + xlab("Mean number of guides")
p <- p + scale_color_viridis(name="target.len")  # viridis color scheme
p <- p + theme_bw(base_size=18) # bw & larger font sizes
p <- p + theme(panel.grid.minor=element_blank(), # leave out minor grid
               aspect.ratio=1)  # square
p <- p + coord_flip()
ggsave("plots/mean-cluster-num-guides-vs-num-seqs-sina.pdf", p, width=8, height=8, useDingbats=FALSE)

# Plot elapsed time vs. number of input seqs
p <- ggplot(summary, aes(x=num.input.seqs, y=elapsed.time.min))
p <- p + geom_point()
p <- p + scale_x_continuous(trans='log10')
p <- p + scale_y_continuous(trans='log10')
p <- p + xlab("Number of sequences") + ylab("Elapsed real time (min)")
p <- p + theme_bw(base_size=18) # bw & larger font sizes
p <- p + theme(panel.grid.minor=element_blank(), # leave out minor grid
               aspect.ratio=1)  # square
ggsave("plots/elapsed-time-vs-num-seqs.pdf", p, width=8, height=8, useDingbats=FALSE)

# Plot memory usage (RSS) vs. number of input seqs
p <- ggplot(summary, aes(x=num.input.seqs, y=rss.mb))
p <- p + geom_point()
p <- p + scale_x_continuous(trans='log10')
p <- p + scale_y_continuous(trans='log10', lim=c(100, 100000))
p <- p + xlab("Number of sequences") + ylab("Memory usage (MB)")
p <- p + theme_bw(base_size=18) # bw & larger font sizes
p <- p + theme(panel.grid.minor=element_blank(), # leave out minor grid
               aspect.ratio=1)  # square
ggsave("plots/memory-vs-num-seqs.pdf", p, width=8, height=8, useDingbats=FALSE)
