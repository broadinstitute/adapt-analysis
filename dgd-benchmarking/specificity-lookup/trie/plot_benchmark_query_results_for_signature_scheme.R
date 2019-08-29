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


# Make a column giving the approach as either 'full' (shard based on 28-mers;
# do not break into two 14-mers and compute signature on full 28-mer) or
# 'split' (shard based on 14-mers; break 28-mer into two 14-mers and compute
# signature from each) or 'noshard' (no sharding; just use a single large,
# regular trie)
# Use 'other' just in case but there should not be any 'other'
# The approach is the prefix of the benchmark value
benchmarks$approach.type <- ifelse(grepl("^full_approach", benchmarks$benchmark),
                                   "full",
                                   ifelse(grepl("^split_approach", benchmarks$benchmark),
                                          "split",
                                          ifelse(grepl("^noshard_approach", benchmarks$benchmark),
                                                 "noshard",
                                                 "other")
                                          )
                                   )

# Remove the name of the approach ('full_approach' or 'split_approach' or
# 'noshard_approach') from the 'benchmark' column
benchmarks$benchmark <- str_replace_all(benchmarks$benchmark, "^(full|split|noshard)_approach_", "")

# Summarize the value across the randomly sampled k-mers for each
# taxonomy and 'experiment' (so there will be a value for each 'experiment'
# for each taxonomy)
benchmarks.summary <- summarySE(benchmarks,
                                measurevar="value",
                                groupvars=c("tax.name", "benchmark",
                                            "approach.type",
                                            "gu.pairing",
                                            "mismatches"))


plot.benchmark.violin <- function(benchmark.name, title, y.lab, no.shard.only, gu.pairing.only,
                                  log10, add1, y.lim, violin.position) {
    # Args:
    #   benchmark.name: name of benchmark to plot ('has.hit', 'num.results',
    #       'nodes_visited', 'runtime')
    #   title: plot title
    #   y.lab: y-axis label
    #   no.shard.only: FALSE/TRUE indicating whether to only show results
    #       from the no-sharding approach
    #   gu.pairing.only: FALSE/TRUE indicating whether to only show results
    #       with GU pairing
    #   log10: FALSE/TRUE indicating whether to plot y-axis on log10 scale
    #   add1: FALSE/TRUE indicating whether to add 1 to all values
    #   y.lim: y-axis limits
    #   violin.position: "identity" for the violin plots to be overlayed
    #       or "dodge" for them to be adjacent to each other

    # Pull out values for the given benchmark
    benchmarks.to.plot <- benchmarks.summary[benchmarks.summary$benchmark == benchmark.name, ]
    benchmarks.to.plot$mismatches <- factor(benchmarks.to.plot$mismatches)

    if (no.shard.only) {
        benchmarks.to.plot <- benchmarks.to.plot[benchmarks.to.plot$approach.type == "noshard", ]
    }
    if (gu.pairing.only) {
        benchmarks.to.plot <- benchmarks.to.plot[benchmarks.to.plot$gu.pairing == "True", ]
    }

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
    p <- p + geom_violin(aes(x=mismatches, y=value, fill=interaction(gu.pairing, approach.type)),
                         alpha=0.5, position=violin.position)
    if (log10) {
        p <- p + scale_y_log10(limits=y.lim)
    } else {
        p <- p + ylim(y.lim)
    }
    p <- p + ggtitle(title)
    p <- p + xlab("Mismatches") + ylab(y.lab)
    p <- p + labs(fill="GU pairing / approach type")
    p <- p + theme_bw()

    return(p)
}


# Make plots with just the no-sharding approach to demonstrate effect
# of GU pairing
p.noshard.has.hit <- plot.benchmark.violin("has_hit", "No sharding: effect of GU pairing on having a hit",
                                           "Fraction of queries with hit",
                                           TRUE, FALSE, FALSE, FALSE, c(0,1), "dodge")
p.noshard.num.results <- plot.benchmark.violin("num_results", "No sharding: effect of GU pairing on number of results",
                                               "1 + Number of results",
                                               TRUE, FALSE, TRUE, TRUE, c(1, 1e4), "dodge")
p.noshard.nodes.visited <- plot.benchmark.violin("nodes_visited", "No sharding: effect of GU pairing on nodes visited",
                                                 "1 + Number of nodes visited",
                                                 TRUE, FALSE, TRUE, TRUE, c(1,1e6), "identity")
p.noshard.runtime <- plot.benchmark.violin("runtime", "No sharding: effect of GU pairing on runtime",
                                           "Runtime (sec)",
                                           TRUE, FALSE, TRUE, FALSE, c(1e-06,10), "identity")

# Make plots with all approaches, but only with GU pairing=TRUE
p.all.has.hit <- plot.benchmark.violin("has_hit", "Sharding approaches: effect on having a hit",
                                       "Fraction of queries with hit",
                                       FALSE, TRUE, FALSE, FALSE, c(0,1), "dodge")
p.all.num.results <- plot.benchmark.violin("num_results", "Sharding approaches: effect on number of results",
                                           "1 + Number of results",
                                           FALSE, TRUE, TRUE, TRUE, c(1, 1e4), "dodge")
p.all.nodes.visited <- plot.benchmark.violin("nodes_visited", "Sharding approaches: effect on number of nodes visited",
                                             "1 + Number of nodes visited",
                                             FALSE, TRUE, TRUE, TRUE, c(1,1e6), "identity")
p.all.runtime <- plot.benchmark.violin("runtime", "Sharding approaches: effect on runtime",
                                       "Runtime (sec)",
                                       FALSE, TRUE, TRUE, FALSE, c(1e-06,10), "identity")

ggsave(OUT.PDF, arrangeGrob(p.noshard.has.hit, p.all.has.hit,
                            p.noshard.num.results, p.all.num.results,
                            p.noshard.nodes.visited, p.all.nodes.visited,
                            p.noshard.runtime, p.all.runtime,
                            ncol=2),
       width=16, height=16, useDingbats=FALSE)

# Remove the empty Rplots.pdf created above
file.remove("Rplots.pdf")
