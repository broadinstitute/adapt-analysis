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
library(ggpubr)
library(reshape2)

# Set threshold for activity deemed highly active
HIGHLY.ACTIVE.THRES <- 2.7198637

summary.fn <- "data/design-results.tsv"
summary <- data.frame(read.table(summary.fn, header=TRUE, sep="\t"))

# Replace '_' in column names with '.'
names(summary) <- gsub("_", ".", names(summary))

# Add some stats
summary$frac.kept <- summary$num.curated.seqs / summary$num.input.seqs

# Add elapsed time and memory in new units
summary$elapsed.time.min <- summary$elapsed.time / 60.0
summary$rss.mb <- summary$rss / 1000.0

# Add column combining taxid and segment
summary$taxseg <- paste0(summary$taxid, "_", summary$segment)

# Make factor grouping num.input.seqs logarithmically
summary$num.input.seqs.group <- ifelse(summary$num.input.seqs < 10, "1-",
                                       ifelse(summary$num.input.seqs < 100, "10-",
                                       ifelse(summary$num.input.seqs < 1000, "100-", ">1000")))
summary$num.input.seqs.group <- factor(summary$num.input.seqs.group, levels=c("1-", "10-", "100-", ">1000"))

# Round mean number of guides, and convert to factor
summary$mean.cluster.num.guides.rounded <- round(summary$mean.cluster.num.guides)
summary[summary$mean.cluster.num.guides.rounded >= 9, ]$mean.cluster.num.guides.rounded <- "ge9"
summary$mean.cluster.num.guides.rounded <- factor(summary$mean.cluster.num.guides.rounded)

# Split by experiment
summary.nonspecific.maxactivity <- summary[summary$experiment == "nonspecific_max-activity", ]
summary.nonspecific.minguides <- summary[summary$experiment == "nonspecific_min-guides", ]
summary.specific.maxactivity <- summary[summary$experiment == "specific_max-activity", ]
summary.specific.minguides <- summary[summary$experiment == "specific_min-guides", ]

# Make data frame comparing expected activity, runtime, etc. (for maximizing activity)
# between nonspecific and specific experiments
shared.taxseg <- intersect(summary.nonspecific.maxactivity$taxseg,
                           summary.specific.maxactivity$taxseg)
summary.compare.specificity <- data.frame(taxseg=shared.taxseg)
summary.compare.specificity$guide.set.expected.activity.nonspecific <- with(summary.nonspecific.maxactivity,
                                                                            mean.cluster.guide.set.expected.activity[match(shared.taxseg, taxseg)])
summary.compare.specificity$guide.set.expected.activity.specific <- with(summary.specific.maxactivity,
                                                                         mean.cluster.guide.set.expected.activity[match(shared.taxseg, taxseg)])
summary.compare.specificity$guide.set.5th.pctile.activity.nonspecific <- with(summary.nonspecific.maxactivity,
                                                                              mean.cluster.guide.set.5th.pctile.activity[match(shared.taxseg, taxseg)])
summary.compare.specificity$guide.set.5th.pctile.activity.specific <- with(summary.specific.maxactivity,
                                                                           mean.cluster.guide.set.5th.pctile.activity[match(shared.taxseg, taxseg)])
summary.compare.specificity$elapsed.time.min.nonspecific <- with(summary.nonspecific.maxactivity,
                                                                 elapsed.time.min[match(shared.taxseg, taxseg)])
summary.compare.specificity$elapsed.time.min.specific <- with(summary.specific.maxactivity,
                                                              elapsed.time.min[match(shared.taxseg, taxseg)])
summary.compare.specificity$rss.mb.nonspecific <- with(summary.nonspecific.maxactivity,
                                                       rss.mb[match(shared.taxseg, taxseg)])
summary.compare.specificity$rss.mb.specific <- with(summary.specific.maxactivity,
                                                    rss.mb[match(shared.taxseg, taxseg)])
summary.compare.specificity$objective.value.nonspecific <- with(summary.nonspecific.maxactivity,
                                                                mean.cluster.objective.value[match(shared.taxseg, taxseg)])
summary.compare.specificity$objective.value.specific <- with(summary.specific.maxactivity,
                                                             mean.cluster.objective.value[match(shared.taxseg, taxseg)])

# Plot histogram of number of clusters per species
# This should be the same for all experiments; use nonspecific.maxactivity
p <- ggplot(summary.nonspecific.maxactivity, aes(num.clusters))
p <- p + geom_histogram(binwidth=1)
p <- p + scale_x_continuous(breaks=seq(1, 15, by=2))
p <- p + xlab("Number of clusters") + ylab("Number of species")
p <- p + theme_pubr()
ggsave("plots/num-clusters.pdf", p, width=4, height=4, useDingbats=FALSE)

# Plot number of clusters vs. log(number of input sequences)
# This should be the same for all experiments; use nonspecific.maxactivity
p <- ggplot(summary.nonspecific.maxactivity, aes(x=num.input.seqs, y=num.clusters))
p <- p + geom_point()
p <- p + scale_x_continuous(trans='log10')
p <- p + scale_y_continuous(breaks=seq(1, 5, by=1))
p <- p + xlab("Number of sequences") + ylab("Number of clusters")
p <- p + theme_pubr()
ggsave("plots/num-clusters-vs-num-seqs.pdf", p, width=4, height=4, useDingbats=FALSE)

# Plot histogram of the fraction of sequences kept after curation
# This should be the same for all experiments; use nonspecific.maxactivity
p <- ggplot(summary.nonspecific.maxactivity, aes(frac.kept))
p <- p + geom_histogram(binwidth=0.05)
p <- p + xlab("Fraction of sequences kept") + ylab("Number of species")
p <- p + theme_pubr()
ggsave("plots/frac-sequences-kept.pdf", p, width=4, height=4, useDingbats=FALSE)

# Plot fraction of sequences kept vs. number of input sequences
# This should be the same for all experiments; use nonspecific.maxactivity
p <- ggplot(summary.nonspecific.maxactivity, aes(x=num.input.seqs, y=frac.kept))
p <- p + geom_point()
p <- p + scale_x_continuous(trans='log10')
p <- p + xlab("Number of sequences") + ylab("Fraction of sequences kept")
p <- p + theme_pubr()
ggsave("plots/frac-sequences-kept-vs-num-seqs.pdf", p, width=4, height=4, useDingbats=FALSE)

# Plot mean genome region (amplicon) length vs. number of input seqs
# (here, means are taken over clusters)
# Use specific.maxactivity experiment
p <- ggplot(summary.specific.maxactivity, aes(x=num.input.seqs, y=mean.cluster.target.len))
p <- p + geom_point(size=1.5, alpha=0.7, shape=16, stroke=0)    # make it easier to see overlapping points
p <- p + scale_x_continuous(trans='log10')
p <- p + xlab("Number of sequences") + ylab("Mean target length")
p <- p + theme_pubr()
ggsave("plots/mean-cluster-target-len-vs-num-seqs.specific_max-activity.pdf", p, width=4, height=4, useDingbats=FALSE)

# Plot mean number of guides vs. number of input seqs
# (here, means are taken over clusters)
# Use specific.maxactivity experiment
p <- ggplot(summary.specific.maxactivity, aes(x=num.input.seqs, y=mean.cluster.num.guides))
p <- p + geom_point(aes(color=mean.cluster.target.len))
p <- p + scale_x_continuous(trans='log10')
p <- p + scale_y_continuous(breaks=seq(1, 20, by=5))
p <- p + xlab("Number of sequences") + ylab("Mean number of guides")
p <- p + scale_color_viridis(name="Length")  # viridis color scheme
p <- p + theme_pubr()
ggsave("plots/mean-cluster-num-guides-vs-num-seqs.specific_max-activity.pdf", p, width=4, height=4, useDingbats=FALSE)

# Plot rounded mean number of guides vs. number of input seqs as sina plot
# Use specific.maxactivity experiment
p <- ggplot(summary.specific.maxactivity, aes(y=num.input.seqs, x=mean.cluster.num.guides.rounded))
p <- p + geom_sina(aes(group=mean.cluster.num.guides.rounded,
                       color=mean.cluster.target.len))
p <- p + scale_y_continuous(trans='log10')
p <- p + ylab("Number of sequences") + xlab("Mean number of guides")
p <- p + scale_color_viridis(name="Length")  # viridis color scheme
p <- p + theme_pubr()
p <- p + coord_flip()
ggsave("plots/mean-cluster-num-guides-vs-num-seqs-sina.specific_max-activity.pdf", p, width=4, height=4, useDingbats=FALSE)

# Plot rounded mean number of guides vs. number of input seqs as sina plot
# Use specific.minguides experiment
p <- ggplot(summary.specific.minguides, aes(y=num.input.seqs, x=mean.cluster.num.guides.rounded))
p <- p + geom_sina(aes(group=mean.cluster.num.guides.rounded,
                       color=mean.cluster.target.len))
p <- p + scale_y_continuous(trans='log10')
p <- p + ylab("Number of sequences") + xlab("Mean number of guides")
p <- p + scale_color_viridis(name="Length")  # viridis color scheme
p <- p + theme_pubr()
p <- p + coord_flip()
ggsave("plots/mean-cluster-num-guides-vs-num-seqs-sina.specific_min-guides.pdf", p, width=4, height=4, useDingbats=FALSE)

# Plot guide set expected activity vs. number of input seqs
# Use specific.maxactivity experiment
p <- ggplot(summary.specific.maxactivity, aes(x=num.input.seqs, y=mean.cluster.guide.set.expected.activity))
p <- p + geom_point(aes(color=mean.cluster.target.len), size=1.5, alpha=0.7, shape=16, stroke=0)
p <- p + scale_x_continuous(trans='log10')
p <- p + xlab("Number of sequences") + ylab("Expected activity of guide set")
p <- p + scale_color_viridis(name="Length")  # viridis color scheme
p <- p + theme_pubr()
ggsave("plots/mean-cluster-guide-set-expected-activity-vs-num-seqs.specific_max-activity.pdf", p, width=4, height=4, useDingbats=FALSE)

# Plot guide set median activity vs. number of input seqs
# Use specific.maxactivity experiment
p <- ggplot(summary.specific.maxactivity, aes(x=num.input.seqs, y=mean.cluster.guide.set.median.activity))
p <- p + geom_point(aes(color=mean.cluster.target.len), size=1.5, alpha=0.7, shape=16, stroke=0)
p <- p + scale_x_continuous(trans='log10')
p <- p + xlab("Number of sequences") + ylab("Median activity of guide set")
p <- p + scale_color_viridis(name="Length")  # viridis color scheme
p <- p + theme_pubr()
ggsave("plots/mean-cluster-guide-set-median-activity-vs-num-seqs.specific_max-activity.pdf", p, width=4, height=4, useDingbats=FALSE)

# Plot guide set 5th percentile activity vs. number of input seqs
# Use specific.maxactivity experiment
p <- ggplot(summary.specific.maxactivity, aes(x=num.input.seqs, y=mean.cluster.guide.set.5th.pctile.activity))
p <- p + geom_point(aes(color=mean.cluster.target.len), size=1.5, alpha=0.7, shape=16, stroke=0)
p <- p + scale_x_continuous(trans='log10')
p <- p + xlab("Number of sequences") + ylab("5th percentile activity of guide set")
p <- p + scale_color_viridis(name="Length")  # viridis color scheme
p <- p + theme_pubr()
ggsave("plots/mean-cluster-guide-set-5th-pctile-activity-vs-num-seqs.specific_max-activity.pdf", p, width=4, height=4, useDingbats=FALSE)

# Plot sina of guide set median and 5th pctile activity, grouped
# by number of input seqs
# Points at 0 overlap, so separate these out to have jitter
# Use specific.maxactivity experiment
summary.specific.maxactivity.melted <- melt(summary.specific.maxactivity, measure.vars=c("mean.cluster.guide.set.median.activity", "mean.cluster.guide.set.5th.pctile.activity"), variable.name="stat", value.name="activity")
summary.specific.maxactivity.melted.nonzero <- summary.specific.maxactivity.melted[summary.specific.maxactivity.melted$activity > 0,]
summary.specific.maxactivity.melted.zero <- summary.specific.maxactivity.melted[summary.specific.maxactivity.melted$activity == 0,]
p <- ggplot(summary.specific.maxactivity.melted.nonzero, aes(x=num.input.seqs.group, y=activity))
p <- p + geom_hline(yintercept=HIGHLY.ACTIVE.THRES, linetype="dashed")
p <- p + geom_sina(aes(group=interaction(num.input.seqs.group, stat),
                       color=stat),
                   size=0.5)
p <- p + geom_point(summary.specific.maxactivity.melted.zero, mapping=aes(x=num.input.seqs.group, y=activity, group=interaction(num.input.seqs.group, stat), color=stat),
                    size=0.5,
                    position=position_jitterdodge())
p <- p + xlab("Number of sequences") + ylab("Activity")
p <- p + theme_pubr()
# Use name="" to avoid legend title
# Use viridis color scheme, but specified manually to avoid the yellow
cols <- c("mean.cluster.guide.set.median.activity"="#643C72",
          "mean.cluster.guide.set.5th.pctile.activity"="#7DBE9C")
labels <- c("mean.cluster.guide.set.median.activity"="Median",
            "mean.cluster.guide.set.5th.pctile.activity"="5th percentile")
p <- p + scale_color_manual(name="",
                            values=cols,
                            labels=labels)
ggsave("plots/mean-cluster-activities-vs-num-seqs-grouped.specific_max-activity.pdf", p, width=4, height=4, useDingbats=FALSE)

# Plot guide set activity for nonspecific vs. specific
p <- ggplot(summary.compare.specificity, aes(x=guide.set.expected.activity.nonspecific, y=guide.set.expected.activity.specific))
p <- p + geom_point(size=1.5, alpha=0.7, shape=16, stroke=0)    # make it easier to see overlapping points
p <- p + geom_abline(slope=1, intercept=0, linetype="dashed")  # diagonal
p <- p + xlim(2, 4) + ylim(2, 4)    # force axes to be the same, but cut out a few outliers
p <- p + xlab("Expected activity of guide set, non-specific") + ylab("Expected activity of guide set, specific")
p <- p + theme_pubr()
ggsave("plots/guide-set-expected-activity-compare-specificity.pdf", p, width=4, height=4, useDingbats=FALSE)

# Plot objective value (primers+guides) for nonspecific vs. specific
p <- ggplot(summary.compare.specificity, aes(x=objective.value.nonspecific, y=objective.value.specific))
p <- p + geom_point(size=1.5, alpha=0.7, shape=16, stroke=0)    # make it easier to see overlapping points
p <- p + geom_abline(slope=1, intercept=0, linetype="dashed")  # diagonal
p <- p + xlim(4, 6) + ylim(4, 6)    # force axes to be the same, but cut out a few outliers
p <- p + xlab("Objective value, non-specific") + ylab("Objective value, specific")
p <- p + theme_pubr()
ggsave("plots/objective-value-compare-specificity.pdf", p, width=4, height=4, useDingbats=FALSE)

# Plot elapsed time vs. number of input seqs
# Use specific.maxactivity experiment
p <- ggplot(summary.specific.maxactivity, aes(x=num.input.seqs, y=elapsed.time.min))
p <- p + geom_point(size=1.5, alpha=0.7, shape=16, stroke=0)    # make it easier to see overlapping points
p <- p + scale_x_continuous(trans='log10')
p <- p + scale_y_continuous(trans='log10')
p <- p + xlab("Number of sequences") + ylab("Elapsed real time (min)")
p <- p + theme_pubr()
ggsave("plots/elapsed-time-vs-num-seqs.specific_max-activity.pdf", p, width=4, height=4, useDingbats=FALSE)

# Plot elapsed time vs. number of input seqs as sina plot
# Use specific.maxactivity experiment
p <- ggplot(summary.specific.maxactivity, aes(x=num.input.seqs.group, y=elapsed.time.min))
p <- p + geom_sina(size=0.5)
p <- p + scale_y_continuous(trans='log10')
p <- p + ylab("Number of sequences") + xlab("Elapsed real time (min)")
p <- p + theme_pubr()
ggsave("plots/elapsed-time-vs-num-seqs-grouped.specific_max-activity.pdf", p, width=4, height=4, useDingbats=FALSE)

ggsave("plots/mean-cluster-num-guides-vs-num-seqs-sina.specific_max-activity.pdf", p, width=4, height=4, useDingbats=FALSE)
# Plot memory usage (RSS) vs. number of input seqs
# Use specific.maxactivity experiment
p <- ggplot(summary.specific.maxactivity, aes(x=num.input.seqs, y=rss.mb))
p <- p + geom_point(size=1.5, alpha=0.7, shape=16, stroke=0)    # make it easier to see overlapping points
p <- p + scale_x_continuous(trans='log10')
p <- p + scale_y_continuous(trans='log10', lim=c(100, 100000))
p <- p + xlab("Number of sequences") + ylab("Memory usage (MB)")
p <- p + theme_pubr()
ggsave("plots/memory-vs-num-seqs.specific_max-activity.pdf", p, width=4, height=4, useDingbats=FALSE)

# Plot elapsed time for nonspecific vs. specific
p <- ggplot(summary.compare.specificity, aes(x=elapsed.time.min.nonspecific, y=elapsed.time.min.specific))
p <- p + geom_point(size=1.5, alpha=0.7, shape=16, stroke=0)    # make it easier to see overlapping points
p <- p + geom_abline(slope=1, intercept=0, linetype="dashed")  # diagonal
p <- p + scale_x_continuous(trans='log10', lim=c(0.1, 10000)) # force axes to be the same
p <- p + scale_y_continuous(trans='log10', lim=c(0.1, 10000))
p <- p + xlab("Elapsed real time, non-specific (min)") + ylab("Elapsed real time, specific (min)")
p <- p + theme_pubr()
ggsave("plots/elapsed-time-compare-specificity.pdf", p, width=4, height=4, useDingbats=FALSE)

# Plot memory usage (RSS) for nonspecific vs. specific
p <- ggplot(summary.compare.specificity, aes(x=rss.mb.nonspecific, y=rss.mb.specific))
p <- p + geom_point(size=1.5, alpha=0.7, shape=16, stroke=0)    # make it easier to see overlapping points
p <- p + geom_abline(slope=1, intercept=0, linetype="dashed")  # diagonal
p <- p + scale_x_continuous(trans='log10', lim=c(100, 100000))  # force axes to be the same
p <- p + scale_y_continuous(trans='log10', lim=c(100, 100000))
p <- p + xlab("Memory usage, non-specific (MB)") + ylab("Memory usage, specific (MB)")
p <- p + theme_pubr()
ggsave("plots/memory-compare-specificity.pdf", p, width=4, height=4, useDingbats=FALSE)
