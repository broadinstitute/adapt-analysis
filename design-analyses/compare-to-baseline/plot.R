#!/usr/bin/env Rscript

# Plot results of comparison against baseline methods.
#
# This uses the results of analyses that were already performed.
#
# Args:
#  1: path to output PDF
#
# By Hayden Metsky <hayden@mit.edu>

library(ggplot2)
library(gridExtra)
library(reshape2)
library(grid)
library(tidyr)


taxonomies <- list.files(path=".", pattern="tax-*")

args <- commandArgs(trailingOnly=TRUE)
out.pdf <- args[1]

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


plot.results.for.taxonomy <- function(taxonomy, real.design.filename,
                                      naive.design.filename, title) {
    # Args:
    #   taxonomy: name of taxonomy being plotted
    #   real.design.filename: name of file containing distribution (one
    #       value per design replicate per window), inside a taxonomy directory
    #       for the "real" designs
    #   naive.design.filename: same as above, for the naive designs
    #   title: title of plot

    # Fill in data frames with a distribution
    real.dist <- data.frame(read.table(
            gzfile(file.path(taxonomy, real.design.filename)), header=TRUE, sep="\t"))
    naive.dist <- data.frame(read.table(
            gzfile(file.path(taxonomy, naive.design.filename)), header=TRUE, sep="\t"))

    # real.dist may be missing windows -- design.py will not output a window
    # if no guides can be constructed for it (e.g., due to missing data) and
    # the window will not show in real.dist if it is missing across all
    # replicates
    # Fill in the missing windows in real.dist
    window.size <- unique(real.dist$window.end - real.dist$window.start)[1]
    real.dist <- complete(real.dist, window.start=full_seq(real.dist$window.start, 1))
    real.dist <- as.data.frame(real.dist)
    real.dist[is.na(real.dist$window.end), ]$window.end <- real.dist[is.na(real.dist$window.end), ]$window.start + window.size

    # Multiply coverage fractions by 100 to obtain percents
    naive.dist$frac.bound.by.consensus <- naive.dist$frac.bound.by.consensus * 100
    naive.dist$frac.bound.by.mode <- naive.dist$frac.bound.by.mode * 100

    # For the real designs, summarize the number of guides in each window
    # across the replicates -- i.e., for each window, find the mean number
    # of guides (and std dev, etc.) across the replicates
    real.dist.summary <- summarySE(real.dist, measurevar="count",
                                   groupvars=c("window.start", "window.end"))
    colnames(real.dist.summary)[colnames(real.dist.summary)=="count"] <- "mean"
    real.dist.summary$approach <- rep("real.guide.count", n=nrow(real.dist.summary))

    # For the naive designs, summarize the coverage obtained by the naive
    # guide in each window across the replicates -- i.e., for each window,
    # find the mean coverage (and std dev, etc.) across the replicates
    # Do this for both the consensus and mode approach
    naive.dist.consensus.summary <- summarySE(naive.dist, measurevar="frac.bound.by.consensus",
                                              groupvars=c("window.start", "window.end"))
    colnames(naive.dist.consensus.summary)[colnames(naive.dist.consensus.summary)=="frac.bound.by.consensus"] <- "mean"
    naive.dist.consensus.summary$approach <- rep("naive.consensus.frac.bound", n=nrow(naive.dist.consensus.summary))
    naive.dist.mode.summary <- summarySE(naive.dist, measurevar="frac.bound.by.mode",
                                         groupvars=c("window.start", "window.end"))
    colnames(naive.dist.mode.summary)[colnames(naive.dist.mode.summary)=="frac.bound.by.mode"] <- "mean"
    naive.dist.mode.summary$approach <- rep("naive.mode.frac.bound", n=nrow(naive.dist.mode.summary))

    # Combine the naive data frames
    naive.dist.summary <- rbind(naive.dist.consensus.summary,
                                naive.dist.mode.summary)
    naive.dist.summary$approach <- factor(naive.dist.summary$approach)

    # Ignore windows where the number of replicates is small, e.g., due
    # to missing data, by making the mean value be NA
    naive.dist.summary[naive.dist.summary$N < 5, ]$mean <- NA
    real.dist.summary[real.dist.summary$N < 5, ]$mean <- NA

    # First produce a plot for the naive designs, showing the fraction
    # of genomes covered at each window
    # Produce plot where x-axis shows the *start* of each window (an
    # alternative would be to compute the middle of each window and use
    # that on the x-axis)
    p1 <- ggplot(naive.dist.summary, aes(x=window.start))
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
    p1.y.min <- quantile(naive.dist.summary$mean - naive.dist.summary$ci,
                         probs=c(0.01), na.rm=TRUE)[1]
    p1.y.max <- quantile(naive.dist.summary$mean + naive.dist.summary$ci,
                         probs=c(0.99), na.rm=TRUE)[1]
    p1 <- p1 + scale_y_continuous(limits=c(p1.y.min, p1.y.max))
    # Add title to plot and axis labels
    p1 <- p1 + ggtitle(title)
    p1 <- p1 + ylab("Coverage against sequences (%)")
    # Leave out usual ggplot2 background and grid lines, but keep border
    # Also leave out the x-axis title/ticks/text
    p1 <- p1 + theme_bw()
    p1 <- p1 + theme(strip.background=element_blank(),
                     panel.grid.major=element_blank(),
                     panel.grid.minor=element_blank(),
                     panel.border=element_rect(colour="black"),
                     plot.title=element_text(size=12),
                     axis.title.x=element_blank(),
                     axis.ticks.x=element_blank(),
                     axis.text.x=element_blank())

    # Second, produce a plot for the "real" design, showing the
    # number of guides required in each window
    # Produce plot where x-axis shows the *start* of each window (an
    # alternative would be to compute the middle of each window and use
    # that on the x-axis)
    p2 <- ggplot(real.dist.summary, aes(x=window.start))
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
    p2.y.min <- quantile(real.dist.summary$mean - real.dist.summary$ci,
                         probs=c(0.01), na.rm=TRUE)[1]
    p2.y.max <- quantile(real.dist.summary$mean + real.dist.summary$ci,
                         probs=c(0.99), na.rm=TRUE)[1]
    p2 <- p2 + scale_y_continuous(limits=c(p2.y.min, p2.y.max))
    # Add axis labels
    p2 <- p2 + xlab("Position") + ylab("Number of guides")
    # Leave out usual ggplot2 background and grid lines, but keep border
    p2 <- p2 + theme_bw()
    p2 <- p2 + theme(strip.background=element_blank(),
                     panel.grid.major=element_blank(),
                     panel.grid.minor=element_blank(),
                     panel.border=element_rect(colour="black"),
                     plot.title=element_blank())

    # Make sure p1 and p2 are on the same x-axis scale (same limit)
    x.min <- min(min(naive.dist.summary$window.start),
                 min(real.dist.summary$window.start))
    x.max <- max(max(naive.dist.summary$window.start),
                 max(real.dist.summary$window.start))
    p1 <- p1 + scale_x_continuous(limits=c(x.min, x.max))
    p2 <- p2 + scale_x_continuous(limits=c(x.min, x.max))

    # Place p1 on top of p2, and draw the plot
    grid.newpage()
    grid.draw(rbind(ggplotGrob(p1), ggplotGrob(p2), size="last"))
}


pdf(out.pdf, width=12, height=8)
for (taxonomy in taxonomies) {
    title <- paste0("Comparison against baseline for ", taxonomy)
    plot.results.for.taxonomy(taxonomy,
                              "real-designs.tsv.gz",
                              "naive-designs.tsv.gz",
                              title)
}
