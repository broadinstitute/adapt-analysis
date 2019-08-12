# Plot results of benchmarking queries with the signature scheme of
# sharding 28-mers into different tries based on computed signature(s).
#
# Args:
#  1: subsampled input size to plot for (i.e., fraction of the total number of
#     k-mers)
#
# By Hayden Metsky <hayden@mit.edu>


require(ggplot2)
require(reshape2)
require(gridExtra)
require(stringr)

args <- commandArgs(trailingOnly=TRUE)
subsampled.size <- args[1]

IN.TABLE <- paste0("out/benchmark-queries-signature-subsample-",
                   subsampled.size,
                   ".tsv.gz")
OUT.PDF <- paste0("out/benchmark-queries-signature-subsample-",
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


# Make a column giving the approach as either 'full' (do not break into
# two 14-mers) or 'split' (break 28-mer into two 14-mers and compute
# signature from each)
# The approach is the prefix of the benchmark value
benchmarks$approach.type <- ifelse(grepl("^full_approach", benchmarks$benchmark), "full", "split")

# Remove the name of the approach ('full_approach' or 'split_approach') from
# the 'benchmark' column
benchmarks$benchmark <- str_replace_all(benchmarks$benchmark, "^(full|split)_approach_", "")

# Summarize the value across the randomly sampled k-mers for each
# taxonomy and 'experiment' (so there will be a value for each 'experiment'
# for each taxonomy)
benchmarks.summary <- summarySE(benchmarks,
                                measurevar="value",
                                groupvars=c("tax.name", "benchmark",
                                            "approach.type",
                                            "mismatches"))


plot.benchmark.violin <- function(benchmark.name, y.lab, log10, add1, y.lim) {
    # Args:
    #   benchmark.name: name of benchmark to plot ('has.hit', 'nodes_visited',
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
    # (colored) for the full vs. split approach
    # Use `position="identity"` so that the full vs. split approaches are
    # shown together (overlayed); otherwise, they would be adjacent (dodge)
    # to each other
    # Note that this assumes gu.pairing is TRUE for all (if not, it will just
    # combine across FALSE/TRUE for gu.pairing)
    p <- p + geom_violin(aes(x=mismatches, y=value, fill=approach.type, color=approach.type),
                         alpha=0.5, position="identity")
    if (log10) {
        p <- p + scale_y_log10(limits=y.lim)
    } else {
        p <- p + ylim(y.lim)
    }
    p <- p + xlab("Mismatches") + ylab(y.lab)
    p <- p + labs(fill="Approach type", color="Approach type")

    return(p)
}


# Make a plot for each benchmark
p.has.hit <- plot.benchmark.violin("has_hit", "Fraction of queries with hit",
                                   FALSE, FALSE, c(0,1))
p.nodes.visited <- plot.benchmark.violin("nodes_visited", "Number of nodes visited",
                                         TRUE, FALSE, c(10,1e6))
p.runtime <- plot.benchmark.violin("runtime", "Runtime (sec)",
                                   TRUE, FALSE, c(1e-06,2))

ggsave(OUT.PDF, arrangeGrob(p.has.hit, p.nodes.visited, p.runtime, ncol=1),
       width=8, height=16, useDingbats=FALSE)

# Remove the empty Rplots.pdf created above
file.remove("Rplots.pdf")
