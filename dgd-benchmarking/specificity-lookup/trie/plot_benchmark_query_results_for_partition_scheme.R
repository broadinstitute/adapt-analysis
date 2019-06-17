# Plot results of benchmarking queries with the partition scheme of
# breaking a 28-mer into two 14-mers.
#
# Args:
#  1: subsampled input size to plot for (i.e., fraction of the total number of
#     k-mers)
#
# By Hayden Metsky <hayden@mit.edu>


require(ggplot2)
require(reshape2)
require(gridExtra)

args <- commandArgs(trailingOnly=TRUE)
subsampled.size <- args[1]

IN.TABLE <- paste0("out/benchmark-queries-partition-14mers-subsample-",
                   subsampled.size,
                   ".tsv.gz")
OUT.PDF <- paste0("out/benchmark-queries-partition-14mers-subsample-",
                  subsampled.size,
                  ".pdf")


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


# Calculate the fraction of all results from the partition scheme that are
# true positives
# The number of reported results from the partition scheme is matches_combined;
# this is actually an upper bound because it sums the number of results from
# the first and second k-mer but there may be overlap between their results --
# a lower bound would be 1/2 this value, so let's uuse 1/2 to get an upperibound
# on the true positive rate
# The number of (correct) results from the 28-mer is matches_from_28_mer
# The true positive rate is
#    [matches_from_28_mer / (0.5*matches_combined)]
# At each point below, lapply will return a list of dataframes; rbind them
# together with do.call(rbind, ...)
tpr.rows <- do.call(rbind, lapply(unique(benchmarks$tax.name),
    function(tax.name) {
        print(paste("Computing TPR for tax.name:", tax.name))
        do.call(rbind, lapply(unique(benchmarks$gu.pairing),
            function(gu.pairing) {
                do.call(rbind, lapply(unique(benchmarks$mismatches),
                    function(mismatches) {
                        vals <- benchmarks[(benchmarks$tax.name == tax.name &
                                            benchmarks$gu.pairing == gu.pairing &
                                            benchmarks$mismatches == mismatches), ]
                        # Compute the true positive rate *for each k-mer* (i.e., row), so that
                        # these will get summarized for each taxonomy and experiment
                        partition.results <- vals[vals$benchmark == "matches_combined", ]$value
                        true.results <- vals[vals$benchmark == "matches_from_28_mer", ]$value
                        tpr <- true.results / (0.5 * partition.results)
                        # nan will result from 0/0; replace these with 0
                        tpr[is.nan(tpr)] <- 0
                        
                        # Make data frame of rows to add
                        num.kmers <- length(tpr)
                        rows.to.add <- data.frame(tax.name=rep(tax.name, num.kmers),
                                                  gu.pairing=rep(gu.pairing, num.kmers),
                                                  mismatches=rep(mismatches, num.kmers),
                                                  benchmark=rep("tpr", num.kmers),
                                                  value=tpr)
                        return(rows.to.add)
                    }
                ))
            }
        ))
    }
))
# Add tpr.rows to benchmarks
benchmarks <- rbind(benchmarks, tpr.rows)

# Summarize the value across the randomly sampled k-mers for each
# taxonomy and 'experiment' (so there will be a value for each 'experiment'
# for each taxonomy)
benchmarks.summary <- summarySE(benchmarks,
                                measurevar="value",
                                groupvars=c("tax.name", "benchmark",
                                            "gu.pairing", "mismatches"))


plot.benchmark.violin <- function(benchmark.name, y.lab, log10, add1, y.lim) {
    # Args:
    #   benchmark.name: name of benchmark to plot ('matches', 'nodes_visited',
    #       'runtime')
    #   y.lab: y-axis label
    #   log10: FALSE/TRUE indicating whether to plot y-axis on log10 scale
    #   add1: FALSE/TRUE indicating whether to add 1 to all values
    #   y.lim: y-axis limits

    # Pull out values for the given benchmark
    benchmarks.to.plot <- benchmarks.summary[benchmarks.summary$benchmark == benchmark.name, ]
    benchmarks.to.plot$mismatches <- factor(benchmarks.to.plot$mismatches)

    if (add1) {
        benchmarks.to.plot$value <- 1 + benchmarks.to.plot$value
    }

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
p.14mers.nodes.visited <- plot.benchmark.violin("nodes_visited", "Number of nodes visited (both 14-mers)", TRUE, FALSE, c(10,1e6))
p.28mers.nodes.visited <- plot.benchmark.violin("nodes_visited_from_28_mer", "Number of nodes visited (28-mer)", TRUE, FALSE, c(10,1e6))
p.14mers.matches <- plot.benchmark.violin("matches_combined", "1 + Number of results (summed across 14-mers)", TRUE, FALSE, c(1,400000))
p.28mers.matches <- plot.benchmark.violin("matches_from_28_mer", "Number of results (28-mer)", FALSE, FALSE, c(0,25))

# It's hard to read a 'tpr' benchmark violin plot because it's so close to 0,
# so show an empirical CDF of it
p.tpr.ecdf <- plot.benchmark.ecdf("tpr", "Fraction of queries with TPR < X", c(0,1))

ggsave(OUT.PDF, arrangeGrob(p.14mers.nodes.visited, p.28mers.nodes.visited, p.14mers.matches, p.28mers.matches, p.tpr.ecdf, ncol=1),
       width=8, height=16, useDingbats=FALSE)

# Remove the empty Rplots.pdf created above
file.remove("Rplots.pdf")
