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
library(ggforce)


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

    # Subset to only keep odd years -- this just makes the data easier
    # to visualize
    dist <- subset(dist, design.year %% 2 == 1 & test.year %% 2 == 1)

    # Only show test.year when it is >= design.year (i.e., do not show
    # test values against years *before* the design was produced, which is
    # not relevant in practice)
    dist <- dist[dist$test.year >= dist$design.year,]

    # For each choice of `design.year` and `test.year` and `sampling` (i.e,
    # bootstrap sample), there are N k-mers (e.g., N=15) with a `frac.hit`
    # value for each
    # Add a column, `kmer.rank`, that gives order of `frac.hit` for each of the
    # N k-mers, within the combination of the 3 variables above -- have the
    # k-mer with the highest `frac.hit` have kmer.rank=1 and the k-mer with
    # the lowest `frac.hit` have kmer.rank=N; when computing rank, negate
    # `frac.hit` to do this
    library(dplyr)
    dist <- dist %>% group_by(design.year, test.year, sampling) %>% mutate(kmer.rank=rank(-frac.hit, ties.method="first"))
    dist <- as.data.frame(dist)
    detach(package:dplyr)

    # Multiply coverage fractions by 100 to obtain percents
    dist$frac.hit <- dist$frac.hit * 100

    # Summarize each coverage value across the replicates -- i.e., for
    # each test.year and design.year and kmer.rank, find
    # the mean coverage (and std dev, etc.) across the samplings (replicates)
    # Effectively, this matches up k-mers across samplings/replicates -- when
    # kmer.rank=1, this will give the mean fraction hit for all the 'best'
    # k-mers from each sampling/replicate and the CI will give a confidence
    # interval across the samplings of frac.hit for the 'best' k-mer, etc.
    dist.summary <- summarySE(dist, measurevar="frac.hit",
                              groupvars=c("design.year", "test.year", "kmer.rank"))

    # Make test.year (facets) be a factor
    dist.summary$test.year.factor <- factor(dist.summary$test.year)

    p <- ggplot(dist.summary,
                aes(color=factor(design.year),
                    x=factor(design.year),
                    y=frac.hit))

    # Plot dots for each combination of design.year and test.year,
    # where each dot represents one of the k-mers for that combination
    # Note this is only incorporating the mean for each k-mer across the
    # bootstrap samples -- it does not show any information about variance
    # across the bootstrap samples; namely, each point is showing
    #  { mean frac.hit across bootstrap samples for the i'th ranked k-mer }
    # Use `method="counts"` to fix a bug and also make the visualization
    # nicer (with `method="density"`, the default), it seems to fail when
    # all values for a design.year/test.year combination are equal (which
    # occurs in one such case, when all values are 0). `scale="count"` also
    # helps the visualization
    p <- p + geom_sina(scale="count",
                       size=1,
                       method="counts",
                       maxwidth=1)

    # Add a horizontal line showing the mean value across k-mers for each
    # combination of design.year and test.year
    p <- p + stat_summary(fun="mean",
                          fun.min="mean",
                          fun.max="mean",
                          size=0.3,   # line thickness
                          geom="crossbar",
                          width=0.8,    # line width
                          show.legend=FALSE)    # do not put crossbar in legend

    # Use viridis color map and label the color legend
    # Use `direction=-1` to reverse colors
    p <- p + scale_color_viridis(discrete=TRUE, name="Design in\nyear",
                                 #direction=-1
                                 )

    # Facet by test.year
    # Allow different scales/space because any test.year may have
    # different number of years in design.year (otherwise, we wind up with
    # lots of empty space in the facets with few different design.year years)
    # Use `switch='x"` to move facet labels to the bottom
    p <- p + facet_grid(. ~ test.year.factor,
                        scales="free_x",
                        space="free_x",
                        switch="x")

    # Make sure the y-axis goes up to 100%; show marks every 10%
    #ci.lower.min <- min(dist.summary$frac.hit - dist.summary$ci)
    #lower.mark <- max(0, round((ci.lower.min - 10)/10.0)*10)
    #p <- p + scale_y_continuous(breaks=seq(lower.mark, 100, 10))
    p <- p + scale_y_continuous(breaks=seq(0, 100, 10))

    # Add axis labels
    # Note xlab() here is actually for the facets, not the x-value in each
    # facet
    p <- p + xlab("Test against year") + ylab("Fraction of sequences with 30-mer (%)")

    # Reformat plot
    p <- p + theme_pubr()

    # Remove x-axis ticks/labels/text (giving design.year) because points
    # are already colored by design.year
    # Also remove the strip background from the facet label
    p <- p + theme(axis.text.x=element_blank(),
                   axis.ticks.x=element_blank(),
                   strip.background=element_blank(),
                   strip.placement="outside",
                   strip.text.x=element_text(size=12))

    return(p)
}

p <- plot.coverage.per.design.per.year(in.tsv)
ggsave(out.pdf, p, width=11.5, height=4, useDingbats=FALSE)
