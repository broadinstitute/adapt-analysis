# Plot results of benchmarking queries.
#
# By Hayden Metsky <hayden@mit.edu>


require(ggplot2)
require(reshape2)
require(gridExtra)

args <- commandArgs(trailingOnly=TRUE)

IN.TABLE <- args[1]
OUT.PDF <- args[2]


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


# Read table and replace '_' in column names with '.'
benchmarks <- read.table(gzfile(IN.TABLE), header=TRUE, sep="\t")
names(benchmarks) <- gsub("_", ".", names(benchmarks))

# Summarize the value across the randomly sampled k-mers for each
# taxonomy and 'experiment' (so there will be a value for each 'experiment'
# for each taxonomy)
benchmarks.summary <- summarySE(benchmarks,
                                measurevar="value",
                                groupvars=c("tax.name", "benchmark",
                                            "gu.pairing", "mismatches"))


plot.benchmark.violin <- function(benchmark.name, y.lab, log10, y.lim) {
    # Args:
    #   benchmark.name: name of benchmark to plot ('matches', 'nodes_visited',
    #       'runtime')
    #   y.lab: y-axis label
    #   log10: FALSE/TRUE indicating whether to plot y-axis on log10 scale
    #   y.lim: y-axis limits

    # Pull out values for the given benchmark
    benchmarks.to.plot <- benchmarks.summary[benchmarks.summary$benchmark == benchmark.name, ]
    benchmarks.to.plot$mismatches <- factor(benchmarks.to.plot$mismatches)

    p <- ggplot(benchmarks.to.plot)

    # Show a violin plot for each value of mismatches, and a separate plot
    # (colored) for with and without G-U pairing
    # Use `position="identity"` so that with and without G-U pairing are
    # shown together (overlayed); otherwise, they would be adjacent (dodge)
    # to each other
    p <- p + geom_violin(aes(x=mismatches, y=value, fill=gu.pairing, color=gu.pairing),
                         alpha=0.5, position="identity")
    if (log10) {
        p <- p + scale_y_log10(limits=y.lim)
    } else {
        p <- p + ylim(y.lim)
    }
    p <- p + xlab("Mismatches") + ylab(y.lab)
    p <- p + labs(fill="G-U pairing", color="G-U pairing")

    return(p)
}


plot.benchmark.ecdf <- function(benchmark.name, y.lab, y.lim) {
    # Args:
    #   benchmark.name: name of benchmark to plot ('matches', 'nodes_visited',
    #       'runtime')
    #   y.lab: y-axis label
    #   y.lim: y-axis limits

    # Pull out values for the given benchmark
    benchmarks.to.plot <- benchmarks.summary[benchmarks.summary$benchmark == benchmark.name, ]
    benchmarks.to.plot$mismatches <- factor(benchmarks.to.plot$mismatches)

    p <- ggplot(benchmarks.to.plot)

    # Show a line for the empirical CDF of the variable; separate line
    # type by G-U pairing and color by mismatches
    p <- p + stat_ecdf(aes(value, linetype=gu.pairing, color=mismatches))

    p <- p + ylim(y.lim)

    p <- p + xlab("X") + ylab(y.lab)
    p <- p + labs(linetype="G-U pairing", color="Mismatches")

    return(p)
}


# Make a plot for each benchmark
p.nodes.visited <- plot.benchmark.violin("nodes_visited", "Number of nodes visited", TRUE, c(10,1e6))
p.runtime <- plot.benchmark.violin("runtime", "Runtime of query (sec)", TRUE, c(1e-5,10))
p.matches <- plot.benchmark.violin("matches", "Number of results", FALSE, c(0,25))


# It's hard to read the 'matches' benchmark plot because almost all
# taxa have 0 mismatches (so the density is heavily concentrated at 0)
# Show the empirical CDF of matches, with a separate line for each
# combination of mismatches and G-U pairing
p.matches.ecdf <- plot.benchmark.ecdf("matches", "Fraction of queries with <= X results", c(0.9, 1.0))

ggsave(OUT.PDF, arrangeGrob(p.nodes.visited, p.runtime, p.matches, p.matches.ecdf, ncol=1),
       width=8, height=12, useDingbats=FALSE)

# Remove the empty Rplots.pdf created above
file.remove("Rplots.pdf")
