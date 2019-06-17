# Plot results of benchmarking queries with the partition scheme of
# breaking a 28-mer into two 14-mers, across input sizes.
#
# Note that plot_benchmark_query_results_for_partition_scheme.R plots
# a true positive rate for the partition scheme, but does not show a plot
# combining results across input sizes.
#
# By Hayden Metsky <hayden@mit.edu>


require(ggplot2)
require(reshape2)
require(gridExtra)

IN.SIZES.STR <- c("0.0001", "0.0002", "0.0004", "0.0008", "0.0016", "0.0032", "0.0064", "0.0128")
IN.TABLE.TEMPLATE <- "out/benchmark-queries-partition-14mers-subsample-__subsampledsize__.tsv.gz"
OUT.PDF <- "out/benchmark-queries-partition-14mers-across-inputs-sizes.pdf"


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


# Construct a table of benchmarks, combining across input sizes
benchmarks <- do.call(rbind, lapply(IN.SIZES.STR,
    function(input.size) {
        # Read table and replace '_' in column names with '.'
        in.table.name <- sub("__subsampledsize__", input.size, IN.TABLE.TEMPLATE)
        benchmarks.for.size <- read.table(gzfile(in.table.name), header=TRUE, sep="\t")
        names(benchmarks.for.size) <- gsub("_", ".", names(benchmarks.for.size))

        # Add column giving the input size
        benchmarks.for.size$input.size <- rep(as.numeric(input.size), nrow(benchmarks.for.size))

        return(benchmarks.for.size)
    }
))


# Summarize the value across the randomly sampled k-mers for each
# input size and 'experiment' (so there will be a value for each 'experiment' and 
# input size)
# Unlike other benchmark plots, here values will be summarized across
# taxonomies
benchmarks.summary <- summarySE(benchmarks,
                                measurevar="value",
                                groupvars=c("input.size", "benchmark",
                                            "gu.pairing", "mismatches"))


# benchmarks is a large table and no longer needed
rm(benchmarks)
gc()


plot.benchmark.for.sizes <- function(benchmark.name, y.lab, log2, add1, y.lim) {
    # Args:
    #   benchmark.name: name of benchmark to plot ('matches', 'nodes_visited',
    #       'runtime')
    #   y.lab: y-axis label
    #   log2: FALSE/TRUE indicating whether to plot y-axis on log2 scale
    #   add1: FALSE/TRUE indicating whether to add 1 to all values
    #   y.lim: y-axis limits
    # Pull out values for the given benchmark
    benchmarks.to.plot <- benchmarks.summary[benchmarks.summary$benchmark == benchmark.name, ]
    benchmarks.to.plot$input.size <- factor(benchmarks.to.plot$input.size)
    benchmarks.to.plot$mismatches <- factor(benchmarks.to.plot$mismatches)

    if (add1) {
        benchmarks.to.plot$value <- 1 + benchmarks.to.plot$value
    }

    p <- ggplot(benchmarks.to.plot)

    # Show a line for each value of mismatches (color) and gu.pairing (line
    # type)
    p <- p + geom_line(aes(x=input.size, y=value,
                           group=interaction(mismatches, gu.pairing),
                           color=mismatches, linetype=gu.pairing))
    if (log2) {
        p <- p + scale_y_continuous(trans='log2', limits=y.lim)
    } else {
        p <- p + ylim(y.lim)
    }
    p <- p + xlab("Fraction of full dataset") + ylab(y.lab)
    p <- p + labs(color="Mismatches", linetype="G-U pairing")

    return(p)
}


# Make a plot for each benchmark
p.14mers.nodes.visited <- plot.benchmark.for.sizes("nodes_visited", "Number of nodes visited (both 14-mers)", TRUE, FALSE, c(10,1e6))
p.28mers.nodes.visited <- plot.benchmark.for.sizes("nodes_visited_from_28_mer", "Number of nodes visited (28-mer)", TRUE, FALSE, c(10,1e6))
p.14mers.matches <- plot.benchmark.for.sizes("matches_combined", "1 + Number of results (summed across 14-mers)", TRUE, TRUE, c(1,400000))
p.28mers.matches <- plot.benchmark.for.sizes("matches_from_28_mer", "1 + Number of results (28-mer)", TRUE, TRUE, c(1,4))

ggsave(OUT.PDF, arrangeGrob(p.14mers.nodes.visited, p.28mers.nodes.visited, p.14mers.matches, p.28mers.matches, ncol=1),
       width=8, height=16, useDingbats=FALSE)

# Remove the empty Rplots.pdf created above
file.remove("Rplots.pdf")
