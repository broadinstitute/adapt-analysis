# Plot comparison of 'greedy' and 'random-greedy' algorithms
# for different choice of parameters.
#
# By Hayden Metsky <hayden@mit.edu>


require(ggplot2)
require(viridis)
require(ggpubr)
require(purrr)

IN.TSV <- "designs.summary.tsv.gz"
OUT.PDF <- "plots/comparison.pdf"


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
designs.summary <- read.table(gzfile(IN.TSV), header=TRUE, sep="\t")
names(designs.summary) <- gsub("_", ".", names(designs.summary))

# Add a column giving the objective value for the guides (not including
# the full target primers and length)
designs.summary$guide.objective.value <- designs.summary$expected.activity - designs.summary$penalty.strength * pmax(0, designs.summary$num.guides - designs.summary$soft.guide.constraint)

# Take the mean across design options for each choice of parameter values
# Do this for guide.objective.value
guide.obj <- summarySE(designs.summary,
                       measurevar="guide.objective.value",
                       groupvars=c("run", "hard.guide.constraint", "soft.guide.constraint", "penalty.strength", "algorithm", "taxonomy"))

# Make facet label names to indicate whether they are
# soft guide constraint (h) or hard guide constraint (H)
soft.guide.constraint.labs <- c("h=1", "h=2", "h=3", "h=4")
names(soft.guide.constraint.labs) <- c(1, 2, 3, 4)
hard.guide.constraint.labs <- c("H=1", "H=2", "H=3", "H=4")
names(hard.guide.constraint.labs) <- c(1, 2, 3, 4)

# Make plots with facets
make.plot <- function(args) {
    t <- args[1]
    ps <- args[2]
    p <- ggboxplot(subset(guide.obj, penalty.strength == ps & taxonomy == t), x="algorithm", y="guide.objective.value", color="algorithm", add="jitter", legend="none") +
        facet_grid(hard.guide.constraint ~ soft.guide.constraint,
                   scales="free_y",
                   labeller=labeller(soft.guide.constraint=soft.guide.constraint.labs, hard.guide.constraint=hard.guide.constraint.labs)) +
        xlab("Algorithm") + ylab("Objective value") +
        theme_pubr() +
        ggtitle(paste0(t, ", lambda=", ps)) +
        rotate_x_text(angle=45) +
        scale_color_viridis(discrete=TRUE)
    return(p)
}

pl <- map(list(c("sars-related-cov", 0.1), c("sars-related-cov", 0.5), c("rhinovirus-a", 0.1), c("rhinovirus-a", 0.5)), make.plot)
g <- ggarrange(plotlist=pl, nrow=2, ncol=2, common.legend=TRUE)
ggsave(OUT.PDF, g, device="pdf", width=16, height=16, useDingbats=FALSE)

# Remove the empty Rplots.pdf created above
file.remove("Rplots.pdf")
