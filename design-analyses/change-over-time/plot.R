#!/usr/bin/env Rscript

# Plot results of analyzing coverage over time.
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


taxonomies <- list.files(path=".", pattern="^tax-*")

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


plot.coverage.per.design.per.year <- function(taxonomy, filename, title) {
    # Args:
    #   taxonomy: name of taxonomy being plotted
    #   filename: name of file containing distribution (one
    #       value per design replicate), inside a taxonomy directory
    #   title: title of plot

    # Fill in a data frame with a distribution, for each design against
    # each year
    dist <- data.frame(read.table(file.path(taxonomy, filename), header=TRUE))

    # Only keep coverage values that are for designs designed using enough
    # sequences AND represent coverage computed against enough sequences
    dist <- dist[(dist$num.seqs.up.to.year >= 10 &
                  dist$num.seqs.against.year >= 10), ]

    # Multiply coverage fractions by 100 to obtain percents
    dist$coverage <- dist$coverage * 100

    # Summarize each coverage value across the replicates -- i.e., for
    # each coverage.against.year and design.up.to.year pair, find
    # the mean coverage (and std dev, etc.) across the replicates
    dist.summary <- summarySE(dist, measurevar="coverage",
                              groupvars=c("coverage.against.year", "design.up.to.year"))

    # Make coverage.against.year (x-axis) be a factor
    dist.summary$coverage.against.year.factor <- factor(dist.summary$coverage.against.year)

    # If there is an 'upto' level, make this the first level
    upto.level <- grep('upto', levels(dist.summary$coverage.against.year.factor), value=TRUE)
    if (length(upto.level) == 1) {
        dist.summary$coverage.against.year.factor <- relevel(dist.summary$coverage.against.year.factor, upto.level)
    }

    p <- ggplot(dist.summary)

    # Plot each point with the confidence interval around each
    p <- p + geom_pointrange(aes(x=coverage.against.year.factor,
                                 y=coverage,
                                 ymin=coverage-ci,
                                 ymax=coverage+ci,
                                 color=factor(design.up.to.year)),
                             size=0.5,
                             fatten=1,
                             position=position_dodge(width=0.5))

    # Be sure to show a label for each value on the x-axis
    #p <- p + scale_x_continuous(breaks=dist.summary$coverage.against.year.factor)

    # Use viridis color map and label the color legend
    p <- p + scale_color_viridis(discrete=TRUE, name="Design in\nyear") 

    # Make sure the y-axis goes up to 100%; show marks every 10%
    ci.lower.min <- min(dist.summary$coverage - dist.summary$ci)
    lower.mark <- max(0, round((ci.lower.min - 10)/10.0)*10)
    p <- p + scale_y_continuous(breaks=seq(lower.mark, 100, 10))

    # Add title to plot and axis labels
    p <- p + ggtitle(title)
    p <- p + xlab("Year X") + ylab("Coverage against sequences from year X (%)")

    # Leave out usual ggplot2 background and grid lines, but keep border
    p <- p + theme_bw()
    p <- p + theme(strip.background=element_blank(),
                   panel.grid.minor=element_blank(),
                   panel.border=element_rect(colour="black"),
                   plot.title=element_text(size=12))

    return(p)
}


taxonomy.plots <- list()
i <- 1
for (taxonomy in taxonomies) {
    title <- paste0("Coverage over time for ", taxonomy)
    p <- plot.coverage.per.design.per.year(taxonomy,
                                           "coverages.distribution.txt",
                                           title)

    # Note that when filling in the taxonomy.plots list, we can do it
    # by explicitly filling in each taxonomy.plots[[i]] but using
    # `append(taxonomy.plots, p)` does not seem to work when p is
    # a plot object
    taxonomy.plots[[i]] <- p
    i <- i + 1
}

g <- grid.arrange(grobs=taxonomy.plots)
ggsave(out.pdf, g, width=8, height=8, useDingbats=FALSE)
