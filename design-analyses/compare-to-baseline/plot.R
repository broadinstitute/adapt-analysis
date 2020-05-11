#!/usr/bin/env Rscript

# Plot results of comparison against baseline methods.
#
# This uses the results of analyses that were already performed.
#
# Args:
#  1: name of taxonomy (must be directory)
#  2: path to output PDF
#
# By Hayden Metsky <hayden@mit.edu>

library(ggplot2)
library(gridExtra)
library(reshape2)
library(grid)
library(tidyr)
library(viridis)
library(ggpubr)


args <- commandArgs(trailingOnly=TRUE)
in.taxonomy <- args[1]
out.pdf <- args[2]

## A helper function from:
##   http://www.cookbook-r.com/Graphs/Plotting_means_and_error_bars_(ggplot2)/#Helper%20functions
## Gives count, mean, standard deviation, standard error of the mean, and
## confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be
##     summarized
##   groupvars: a vector containing names of columns that contain grouping
##     variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is
##     95%)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
    library(plyr)

    # New version of length which can handle NA's: if na.rm==T, don't count
    # them
    length2 <- function (x, na.rm=FALSE) {
        if (na.rm) sum(!is.na(x))
        else       length(x)
    }

    # This does the summary. For each group's data frame, return a vector with
    # N, mean, and sd
    datac <- ddply(data, groupvars, .drop=.drop,
      .fun = function(xx, col) {
        c(N    = length2(xx[[col]], na.rm=na.rm),
          mean = mean   (xx[[col]], na.rm=na.rm),
          sd   = sd     (xx[[col]], na.rm=na.rm)
        )
      },
      measurevar
    )

    # Rename the "mean" column    
    datac <- rename(datac, c("mean" = measurevar))

    datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean

    # Confidence interval multiplier for standard error
    # Calculate t-statistic for confidence interval: 
    # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
    ciMult <- qt(conf.interval/2 + .5, datac$N-1)
    datac$ci <- datac$se * ciMult

    return(datac)
}


plot.results.for.taxonomy <- function(taxonomy) {
    # Args:
    #   taxonomy: name of taxonomy being plotted

    # Fill in data frames with a distribution
    real.max.dist <- data.frame(read.table(
            gzfile(file.path(taxonomy, "real-designs.max-activity.tsv.gz")),
            header=TRUE,
            sep="\t", na.strings=c("NA", "None")))
    real.min.dist <- data.frame(read.table(
            gzfile(file.path(taxonomy, "real-designs.min-guides.tsv.gz")),
            header=TRUE,
            sep="\t", na.strings=c("NA", "None")))
    naive.dist <- data.frame(read.table(
            gzfile(file.path(taxonomy, "naive-designs.tsv.gz")), header=TRUE,
            sep="\t", na.strings=c("NA", "None")))

    # If the guide sequence is NA, then interpret the count or coverage as NA too
    # (it may be 0)
    real.min.dist$count[which(is.na(real.min.dist$target.sequences))] <- NA
    real.min.dist$total.frac.bound[which(is.na(real.min.dist$target.sequences))] <- NA
    real.min.dist$score[which(is.na(real.min.dist$target.sequences))] <- NA
    real.max.dist$count[which(is.na(real.max.dist$target.sequences))] <- NA
    real.max.dist$total.frac.bound[which(is.na(real.max.dist$target.sequences))] <- NA
    naive.dist$frac.bound.by.consensus[which(is.na(naive.dist$target.sequence.by.consensus))] <- NA
    naive.dist$frac.bound.by.mode[which(is.na(naive.dist$target.sequence.by.mode))] <- NA

    # real.*.dist may be missing windows -- design.py will not output a window
    # if no guides can be constructed for it (e.g., due to missing data) and
    # the window will not show in real.*.dist if it is missing across all
    # replicates
    # Fill in the missing windows in real.*.dist
    window.size <- unique(real.max.dist$window.end - real.max.dist$window.start)[1]
    real.max.dist <- complete(real.max.dist, window.start=full_seq(real.max.dist$window.start, 1))
    real.max.dist <- as.data.frame(real.max.dist)
    real.max.dist$window.end[which(is.na(real.max.dist$window.end))] <- real.max.dist$window.start[which(is.na(real.max.dist$window.end))] + window.size
    real.min.dist <- complete(real.min.dist, window.start=full_seq(real.min.dist$window.start, 1))
    real.min.dist <- as.data.frame(real.min.dist)
    real.min.dist$window.end[which(is.na(real.min.dist$window.end))] <- real.min.dist$window.start[which(is.na(real.min.dist$window.end))] + window.size

    # Multiply coverage fractions by 100 to obtain percents
    naive.dist$frac.bound.by.consensus <- naive.dist$frac.bound.by.consensus * 100
    naive.dist$frac.bound.by.mode <- naive.dist$frac.bound.by.mode * 100
    real.max.dist$frac.bound <- real.max.dist$total.frac.bound * 100

    # For the real.min designs, summarize the number of guides in each window
    # across the replicates -- i.e., for each window, find the mean number
    # of guides (and std dev, etc.) across the replicates
    real.min.dist.summary <- summarySE(real.min.dist, measurevar="count",
                                   groupvars=c("window.start", "window.end"))
    colnames(real.min.dist.summary)[colnames(real.min.dist.summary)=="count"] <- "mean"
    real.min.dist.summary$approach <- "real.min.guide.count"

    # For the real.max designs, summarize the coverage obtained by the
    # guides in each window and choice of hard guide constraint (hgc) across
    # the replicates -- i.e., for each window and hgc,
    # find the mean coverage (and std dev, etc.) across the replicates
    real.max.dist.summary <- summarySE(real.max.dist, measurevar="frac.bound",
                                              groupvars=c("window.start", "window.end", "hgc"))
    colnames(real.max.dist.summary)[colnames(real.max.dist.summary)=="frac.bound"] <- "mean"
    real.max.dist.summary$approach <- "real.max.activity.frac.bound"

    # For the naive designs, summarize the coverage obtained by the naive
    # guide in each window across the replicates -- i.e., for each window,
    # find the mean coverage (and std dev, etc.) across the replicates
    # Do this for both the consensus and mode approach
    naive.dist.consensus.summary <- summarySE(naive.dist, measurevar="frac.bound.by.consensus",
                                              groupvars=c("window.start", "window.end"))
    colnames(naive.dist.consensus.summary)[colnames(naive.dist.consensus.summary)=="frac.bound.by.consensus"] <- "mean"
    naive.dist.consensus.summary$approach <- "naive.consensus.frac.bound"
    naive.dist.mode.summary <- summarySE(naive.dist, measurevar="frac.bound.by.mode",
                                         groupvars=c("window.start", "window.end"))
    colnames(naive.dist.mode.summary)[colnames(naive.dist.mode.summary)=="frac.bound.by.mode"] <- "mean"
    naive.dist.mode.summary$approach <- "naive.mode.frac.bound"

    # Combine the naive data frames
    naive.dist.summary <- rbind(naive.dist.consensus.summary,
                                naive.dist.mode.summary)
    naive.dist.summary$approach <- factor(naive.dist.summary$approach)

    # For the real design max-activity data, only show where the number of
    # guides (hgc -- i.e., hard guide constraint) is 1 or 2
    real.max.dist.summary <- real.max.dist.summary[real.max.dist.summary$hgc %in% c(1,2),]

    # Combine naive.dist.summary with real.max.dist.summary -- both report
    # coverages (frac.bound)
    real.max.dist.summary$approach <- paste0(real.max.dist.summary$approach, ".hgc-",
                                             real.max.dist.summary$hgc)
    real.max.dist.summary <- subset(real.max.dist.summary, select=-c(hgc))
    frac.bounds.summary <- rbind(naive.dist.summary, real.max.dist.summary)

    # Ignore windows where the number of replicates is small (e.g., due
    # to missing data or too many gaps) by making the mean value be NA
    frac.bounds.summary$mean[which(frac.bounds.summary$N < 5)] <- NA
    real.min.dist.summary$mean[which(real.min.dist.summary$N < 5)] <- NA

    # First produce a plot for the naive and real.max-activity designs, showing the fraction
    # of genomes covered at each window
    # Produce plot where x-axis shows the *start* of each window (an
    # alternative would be to compute the middle of each window and use
    # that on the x-axis)
    p1 <- ggplot(frac.bounds.summary, aes(x=window.start))
    # Plot mean values as a line
    p1 <- p1 + geom_line(aes(y=mean, color=approach), size=1.5)
    # Can use geom_errorbar(..) to show error bars at each plotted
    # x-value; alternatively, geom_ribbon(..) to show a continuous
    # interval (i.e., confidence band) around each line. Note that
    # this is a 95% pointwise confidence band, NOT a simultaneous
    # confidence band.
    p1 <- p1 + geom_ribbon(aes(ymin=mean-ci,
                               ymax=mean+ci,
                               fill=approach), alpha=0.2)
    # Manually set the y-axis limits to avoid outliers in the confidence
    # intervals (the ribbon); these outliers may not show on the plot
    # (use the outlier definition of any point <(Qlo - 1.5*IQR) or
    # >(Qhi + 1.5*IQR) where IQR=Qhi-Qlo and Qlo=10% pctile and Qhi=90% pctile,
    # and the Qlo/Qhi are computed across the genome based on the mean +/- ci
    # value at each position)
    p1.y.Qlo <- quantile(frac.bounds.summary$mean - frac.bounds.summary$ci,
                         probs=c(0.1), na.rm=TRUE)[1]
    p1.y.Qhi <- quantile(frac.bounds.summary$mean + frac.bounds.summary$ci,
                         probs=c(0.9), na.rm=TRUE)[1]
    p1.y.min <- max(min(frac.bounds.summary$mean - frac.bounds.summary$ci),
                    p1.y.Qlo - 1.5*(p1.y.Qhi - p1.y.Qlo))
    p1.y.max <- min(max(frac.bounds.summary$mean + frac.bounds.summary$ci),
                    p1.y.Qhi + 1.5*(p1.y.Qhi - p1.y.Qlo))
    p1 <- p1 + scale_y_continuous(limits=c(p1.y.min, p1.y.max))
    # Add title to plot and axis labels
    p1 <- p1 + ggtitle(taxonomy)
    p1 <- p1 + xlab("Position") + ylab("Coverage against sequences (%)")
    p1 <- p1 + theme_pubr()

    # Second, produce a plot for the real.min-guides design, showing the
    # number of guides required in each window
    # Produce plot where x-axis shows the *start* of each window (an
    # alternative would be to compute the middle of each window and use
    # that on the x-axis)
    p2 <- ggplot(real.min.dist.summary, aes(x=window.start))
    # Plot mean values as points
    p2 <- p2 + geom_line(aes(y=mean, color=approach), size=1.5)
    # Can use geom_errorbar(..) to show error bars at each plotted
    # x-value; alternatively, geom_ribbon(..) to show a continuous
    # interval (i.e., confidence band) around each line. Note that
    # this is a 95% pointwise confidence band, NOT a simultaneous
    # confidence band.
    p2 <- p2 + geom_ribbon(aes(ymin=mean-ci,
                               ymax=mean+ci,
                               fill=approach), alpha=0.2)
    # Manually set the y-axis limits to avoid outliers in the confidence
    # intervals (the ribbon); these outliers may not show on the plot
    p2.y.Qlo <- quantile(real.min.dist.summary$mean - real.min.dist.summary$ci,
                         probs=c(0.1), na.rm=TRUE)[1]
    p2.y.Qhi <- quantile(real.min.dist.summary$mean + real.min.dist.summary$ci,
                         probs=c(0.9), na.rm=TRUE)[1]
    p2.y.min <- max(min(real.min.dist.summary$mean - real.min.dist.summary$ci),
                    p2.y.Qlo - 1.5*(p2.y.Qhi - p2.y.Qlo))
    p2.y.max <- min(max(real.min.dist.summary$mean + real.min.dist.summary$ci),
                    p2.y.Qhi + 1.5*(p2.y.Qhi - p2.y.Qlo))
    p2.y.min <- min(0, p2.y.min)    # make sure to show y=0 guides
    p2 <- p2 + scale_y_continuous(limits=c(p2.y.min, p2.y.max))
    # Add axis labels
    p2 <- p2 + xlab("Position") + ylab("Number of guides")
    p2 <- p2 + theme_pubr()

    # Make sure p1 and p2 are on the same x-axis scale (same limit)
    x.min <- min(min(frac.bounds.summary$window.start),
                 min(real.min.dist.summary$window.start))
    x.max <- max(max(frac.bounds.summary$window.start),
                 max(real.min.dist.summary$window.start))
    p1 <- p1 + scale_x_continuous(limits=c(x.min, x.max))
    p2 <- p2 + scale_x_continuous(limits=c(x.min, x.max))

    # Use viridis color map for both plots
    p1 <- p1 + scale_color_viridis(discrete=TRUE)
    p1 <- p1 + scale_fill_viridis(discrete=TRUE)
    p2 <- p2 + scale_color_viridis(discrete=TRUE)
    p2 <- p2 + scale_fill_viridis(discrete=TRUE)

    # Place p1 on top of p2
    g <- grid.arrange(p1, p2, ncol=1)
    return(g)
}


g <- plot.results.for.taxonomy(in.taxonomy)
ggsave(out.pdf, g, width=12, height=8, useDingbats=FALSE)
