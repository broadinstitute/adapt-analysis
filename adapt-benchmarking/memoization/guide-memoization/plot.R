# Plot elapsed runtime vs. window during search with and without memoization.
#
# By Hayden Metsky <hayden@mit.edu>


require(ggplot2)
require(gridExtra)
require(reshape2)
require(viridis)
require(ggpubr)

args <- commandArgs(trailingOnly=TRUE)
IN.TSV <- args[1]
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


# Read TSV of results and replace '_' in column names with '.'
runtime <- read.table(gzfile(IN.TSV), header=TRUE, sep="\t")
names(runtime) <- gsub("_", ".", names(runtime))

# The time is in milliseconds; convert it to minutes
runtime$time.from.start <- runtime$time.from.start / 1000.0 / 60.0

# To cut down on the number of plots to point, only use every 10th window
# This will also leave out the 1st window, which helps when plotting on
# the y-axis (since the time is, by definition, 0 at the first window)
runtime <- subset(runtime, window.num %% 10 == 0)

# Remove window.start and window.end, which we don't use
runtime$window.start <- NULL
runtime$window.end <- NULL

# To determine how many windows to plot, find the maximum window.num
# for each group of (run, memoize) -- and then take the minimum
# across all of these
require(dplyr)
max.windows.by.group <- runtime %>%
    group_by(run, memoize) %>% summarise(window.num=max(window.num))
max.window.num <- min(max.windows.by.group$window.num)
runtime <- subset(runtime, window.num <= max.window.num)

# Summarize over the folds (splits) of the training data by grouping
# based on all the other variables
# Summarize over different random runs by grouping on other variables
runtime.summarized <- summarySE(runtime,
                                measurevar="time.from.start",
                                groupvars=c("memoize", "window.num"))

# Rename memoize factor levels
runtime.summarized$memoize <- factor(runtime.summarized$memoize, levels=c("no", "yes"))
runtime.summarized$memoize <- recode(runtime.summarized$memoize,
                                     "no"="Without memoization",
                                     "yes"="With memoization")

# Calculate ribbon hi/lo, but keep lo > .1 so it doesn't go negative,
# which would cause problems on a log scale
runtime.summarized$time.from.start.hi <- runtime.summarized$time.from.start + runtime.summarized$ci
runtime.summarized$time.from.start.lo <- pmax(0.1, runtime.summarized$time.from.start - runtime.summarized$ci)

# Make plot
memoize.colors <- c("With memoization"="#6CC199", "Without memoization"="#6C3977")
p <- ggplot(runtime.summarized, aes(x=window.num, y=time.from.start))
p <- p + geom_line(aes(color=memoize))
p <- p + geom_ribbon(aes(ymin=time.from.start.lo, ymax=time.from.start.hi, fill=memoize), alpha=0.3)
p <- p + scale_y_log10()    # log axis
p <- p + xlab("Number of amplicons searched") + ylab("Elapsed real time (min)")
p <- p + scale_fill_manual(values=memoize.colors)
p <- p + scale_color_manual(values=memoize.colors)
p <- p + labs(color="", fill="")    # no legend title
p <- p + theme_pubr(legend="top")

ggsave(OUT.PDF, p, width=4, height=4, useDingbats=FALSE)

