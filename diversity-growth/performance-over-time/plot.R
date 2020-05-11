#!/usr/bin/env Rscript

# Plot results of coverage over time.
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
library(viridis)
library(ggpubr)


args <- commandArgs(trailingOnly=TRUE)
in.tsv <- args[1]
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


plot.coverage.per.design.per.year <- function(data.path) {
    # Args:
    #   data.path: path to data TSV

    # Fill in a data frame with a distribution, for each design against
    # each year
    dist <- data.frame(read.table(data.path, sep="\t", header=TRUE))

    # Replace '_' in column names with '.'
    names(dist) <- gsub("_", ".", names(dist))

    # Only show test.year when it is >= design.year (i.e., do not show
    # test values against years *before* the design was produced, which is
    # not relevant in practice)
    dist <- dist[dist$test.year >= dist$design.year,]

    # Multiply coverage fractions by 100 to obtain percents
    dist$frac.hit <- dist$frac.hit * 100

    # Summarize each coverage value across the replicates -- i.e., for
    # each test.year and design.year pair, find
    # the mean coverage (and std dev, etc.) across the replicates
    dist.summary <- summarySE(dist, measurevar="frac.hit",
                              groupvars=c("design.year", "test.year"))

    # Make test.year (x-axis) be a factor
    dist.summary$test.year.factor <- factor(dist.summary$test.year)

    p <- ggplot(dist.summary)

    # Plot each point with the confidence interval around each
    p <- p + geom_pointrange(aes(x=test.year.factor,
                                 y=frac.hit,
                                 ymin=frac.hit-ci,
                                 ymax=frac.hit+ci,
                                 color=factor(design.year)),
                             size=0.5,
                             fatten=1,
                             position=position_dodge(width=0.4))

    # Use viridis color map and label the color legend
    p <- p + scale_color_viridis(discrete=TRUE, name="Design in\nyear")

    # Make sure the y-axis goes up to 100%; show marks every 10%
    ci.lower.min <- min(dist.summary$frac.hit - dist.summary$ci)
    lower.mark <- max(0, round((ci.lower.min - 10)/10.0)*10)
    p <- p + scale_y_continuous(breaks=seq(lower.mark, 100, 10))

    # Add title to plot and axis labels
    p <- p + xlab("Year X") + ylab("Average fraction of sequences from year X detected by most conserved 30-mers (%)")

    # Reformat plot
    p <- p + theme_pubr()

    return(p)
}

p <- plot.coverage.per.design.per.year(in.tsv)
ggsave(out.pdf, p, width=16, height=16, useDingbats=FALSE)
