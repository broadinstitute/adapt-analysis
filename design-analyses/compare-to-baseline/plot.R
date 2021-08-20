#!/usr/bin/env Rscript

# Plot results of comparison against baseline methods.
#
# This uses the results of analyses that were already performed.
#
# Note that there may be empty (not plotted) regions of the genome
# in the minimize-guides plot (bottom plot -- i.e, p2). This represents
# regions of the genome for which no number of guides could achieve
# the coverage threshold, for example, due to missing data or gaps.
#
# Args:
#  1: name of taxonomy (must be directory)
#  2: path to output directory in which to place PDFs
#
# By Hayden Metsky <hayden@mit.edu>

library(ggplot2)
library(gridExtra)
library(reshape2)
library(grid)
library(tidyr)
library(viridis)
library(ggpubr)
library(plyr)


args <- commandArgs(trailingOnly=TRUE)
in.taxonomy <- args[1]
taxonomy.name <- args[2]
out.dir <- args[3]

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


plot.results.for.taxonomy <- function(taxonomy) {
    # Args:
    #   taxonomy: name of taxonomy being plotted

    # Fill in data frames with a distribution
    real.max.dist <- data.frame(read.table(
            gzfile(file.path(taxonomy, "real-designs.max-activity.tsv.gz")),
            header=TRUE,
            sep="\t", na.strings=c("NA", "None")))
    real.min.dist <- data.frame(read.table(
            gzfile(file.path(taxonomy, "real-designs.min-guides.tsv.gz")),
            header=TRUE,
            sep="\t", na.strings=c("NA", "None")))
    naive.dist <- data.frame(read.table(
            gzfile(file.path(taxonomy, "naive-designs.tsv.gz")), header=TRUE,
            sep="\t", na.strings=c("NA", "None")))

    # If the guide sequence is NA, then interpret the count or coverage as NA too
    # (it may be 0)
    real.min.dist$count[which(is.na(real.min.dist$target.sequences))] <- NA
    real.min.dist$total.frac.bound[which(is.na(real.min.dist$target.sequences))] <- NA
    real.min.dist$score[which(is.na(real.min.dist$target.sequences))] <- NA
    real.max.dist$count[which(is.na(real.max.dist$target.sequences))] <- NA
    real.max.dist$total.frac.bound[which(is.na(real.max.dist$target.sequences))] <- NA
    naive.dist$frac.bound[which(is.na(naive.dist$target.sequence))] <- NA

    # real.*.dist may be missing windows -- design.py will not output a window
    # if no guides can be constructed for it (e.g., due to missing data) and
    # the window will not show in real.*.dist if it is missing across all
    # replicates
    # Fill in the missing windows in real.*.dist
    window.size <- unique(real.max.dist$window.end - real.max.dist$window.start)[1]
    #real.max.dist <- complete(real.max.dist, window.start=full_seq(real.max.dist$window.start, 1))
    #real.max.dist <- as.data.frame(real.max.dist)
    #real.max.dist$window.end[which(is.na(real.max.dist$window.end))] <- real.max.dist$window.start[which(is.na(real.max.dist$window.end))] + window.size
    #real.min.dist <- complete(real.min.dist, window.start=full_seq(real.min.dist$window.start, 1))
    #real.min.dist <- as.data.frame(real.min.dist)
    #real.min.dist$window.end[which(is.na(real.min.dist$window.end))] <- real.min.dist$window.start[which(is.na(real.min.dist$window.end))] + window.size

    # Multiply coverage fractions by 100 to obtain percents
    naive.dist$frac.bound <- naive.dist$frac.bound * 100
    real.max.dist$frac.bound <- real.max.dist$total.frac.bound * 100

    # For the real.min designs, summarize the number of guides in each window
    # and choice of guide coverage (gp) across the replicates -- i.e., for
    # each window, find the mean number of guides (and std dev, etc.) across
    # the replicates
    real.min.dist.summary <- summarySE(real.min.dist, measurevar="count",
                                   groupvars=c("window.start", "window.end", "gp"))
    colnames(real.min.dist.summary)[colnames(real.min.dist.summary)=="count"] <- "mean"
    real.min.dist.summary$approach <- "real.min.guide.count"

    # For the real.max designs, summarize the coverage obtained by the
    # guides in each window and choice of hard guide constraint (hgc) across
    # the replicates -- i.e., for each window and hgc,
    # find the mean coverage (and std dev, etc.) across the replicates
    real.max.dist.summary <- summarySE(real.max.dist, measurevar="frac.bound",
                                              groupvars=c("window.start", "window.end", "hgc"))
    colnames(real.max.dist.summary)[colnames(real.max.dist.summary)=="frac.bound"] <- "mean"
    real.max.dist.summary$approach <- "real.max.activity.frac.bound"

    # For the naive designs, summarize the coverage obtained by the naive
    # guide in each window across the replicates -- i.e., for each window,
    # find the mean coverage (and std dev, etc.) across the replicates
    # Group by window, approach (consensus or mode), and number of probes
    naive.dist.summary <- summarySE(naive.dist, measurevar="frac.bound",
                                    groupvars=c("window.start", "window.end", "approach", "num.probes"))
    colnames(naive.dist.summary)[colnames(naive.dist.summary)=="frac.bound"] <- "mean"
    naive.dist.summary$approach <- as.character(naive.dist.summary$approach)
    naive.dist.summary$approach[naive.dist.summary$approach=="consensus"] <- "naive.consensus.frac.bound"
    naive.dist.summary$approach[naive.dist.summary$approach=="mode"] <- paste0("naive.mode.frac.bound", ".upto-",
                                                                               naive.dist.summary$num.probes[naive.dist.summary$approach=="mode"])
    naive.dist.summary$approach <- factor(naive.dist.summary$approach)

    # For the real design max-activity data, only show where the number of
    # guides (hgc -- i.e., hard guide constraint) is 1 or 2 or 3
    real.max.dist.summary <- real.max.dist.summary[real.max.dist.summary$hgc %in% c(1,2,3),]

    # Combine naive.dist.summary with real.max.dist.summary -- both report
    # coverages (frac.bound)
    real.max.dist.summary$approach <- paste0(real.max.dist.summary$approach, ".hgc-",
                                             real.max.dist.summary$hgc)
    real.max.dist.summary <- subset(real.max.dist.summary, select=-c(hgc))
    naive.dist.summary <- subset(naive.dist.summary, select=-c(num.probes))
    frac.bounds.summary <- rbind(naive.dist.summary, real.max.dist.summary)

    # In real.min.dist.summary, make guide coverage (gp) be a percent and factor
    real.min.dist.summary$gp.display <- paste0(real.min.dist.summary$gp * 100, "%")
    real.min.dist.summary$gp.display <- factor(real.min.dist.summary$gp.display)

    # Ignore windows where the number of replicates is small (e.g., due
    # to missing data or too many gaps) by making the mean value be NA
    frac.bounds.summary$mean[which(frac.bounds.summary$N < 5)] <- NA
    real.min.dist.summary$mean[which(real.min.dist.summary$N < 5)] <- NA

    # Remove windows too close to the start or end of the alignment (within
    # 10 nt) -- because genomes are of different lengths, there may be large
    # drops in coverage at the ends
    aln.end <- max(frac.bounds.summary$window.end)
    len.tol <- 10
    frac.bounds.summary <- frac.bounds.summary[frac.bounds.summary$window.start >= 0 + len.tol,]
    frac.bounds.summary <- frac.bounds.summary[frac.bounds.summary$window.end < aln.end - len.tol,]
    real.min.dist.summary <- real.min.dist.summary[real.min.dist.summary$window.start >= 0 + len.tol,]
    real.min.dist.summary <- real.min.dist.summary[real.min.dist.summary$window.end < aln.end - len.tol,]

    # Use a window midpoint as position
    frac.bounds.summary$aln.pos <- frac.bounds.summary$window.start + (frac.bounds.summary$window.end - frac.bounds.summary$window.start)/2
    real.min.dist.summary$aln.pos <- real.min.dist.summary$window.start + (real.min.dist.summary$window.end - real.min.dist.summary$window.start)/2

    # Map the positions from alignment-space to a position in a reference
    # genome
    ref.pos.map <- as.data.frame(read.table(file.path(taxonomy, "input-alns/ref-pos-map.tsv"), sep="\t", header=TRUE))
    frac.bounds.summary <- merge(frac.bounds.summary, ref.pos.map, by="aln.pos")
    real.min.dist.summary <- merge(real.min.dist.summary, ref.pos.map, by="aln.pos")

    # Split frac.bounds.summary into separate data frames, for different
    # types of plots:
    #   1) frac.bounds.summary will have the consensus approach and the mode
    #      approach with only 1 probe
    #   2) frac.bounds.summary.all.mode will *not* have the consensus approach
    #      but will have the mode approach with all numbers of probes
    # Both will keep the real.max.dist.summary (i.e., ADAPT) rows
    frac.bounds.summary.all.mode <- frac.bounds.summary[frac.bounds.summary$approach != "naive.consensus.frac.bound",]
    frac.bounds.summary <- frac.bounds.summary[!(frac.bounds.summary$approach %in% paste0("naive.mode.frac.bound.upto-", c(2:10))),]

    # Rename approaches to display
    frac.bounds.summary$approach.display <- mapvalues(frac.bounds.summary$approach,
        from=c("naive.consensus.frac.bound", "naive.mode.frac.bound.upto-1", "real.max.activity.frac.bound.hgc-1", "real.max.activity.frac.bound.hgc-2", "real.max.activity.frac.bound.hgc-3"),
        to=c("Consensus", "Mode", "ADAPT, 1 probe", "ADAPT, 2 probes", "ADAPT, 3 probes"))
    frac.bounds.summary.all.mode$approach.display <- mapvalues(frac.bounds.summary.all.mode$approach,
        from=c(paste0("naive.mode.frac.bound.upto-", c(1:10)), "real.max.activity.frac.bound.hgc-1", "real.max.activity.frac.bound.hgc-2", "real.max.activity.frac.bound.hgc-3"),
        to=c(paste0("Mode (", c(1:10), ")"), "ADAPT, 1 probe", "ADAPT, 2 probes", "ADAPT, 3 probes"))
    real.min.dist.summary$approach.display <- mapvalues(real.min.dist.summary$approach,
        from=c("real.min.guide.count"),
        to=c("Minimize probes"))

    # First produce a plot for the naive and real.max-activity designs, showing the fraction
    # of genomes covered at each window
    # Produce plot where x-axis shows the reference position of the middle of
    # the window
    p1 <- ggplot(frac.bounds.summary, aes(x=ref.pos))
    # Plot mean values as a line
    p1 <- p1 + geom_line(aes(y=mean, color=approach.display), size=1.0)
    # Can use geom_errorbar(..) to show error bars at each plotted
    # x-value; alternatively, geom_ribbon(..) to show a continuous
    # interval (i.e., confidence band) around each line. Note that
    # this is a 95% pointwise confidence band, NOT a simultaneous
    # confidence band.
    p1 <- p1 + geom_ribbon(aes(ymin=mean-ci,
                               ymax=mean+ci,
                               fill=approach.display), alpha=0.2)
    # Manually set the y-axis limits to avoid outliers in the confidence
    # intervals (the ribbon); these outliers may not show on the plot
    # (use the outlier definition of any point <(Qlo - 1.5*IQR) or
    # >(Qhi + 1.5*IQR) where IQR=Qhi-Qlo and Qlo=10% pctile and Qhi=90% pctile,
    # and the Qlo/Qhi are computed across the genome based on the mean +/- ci
    # value at each position)
    p1.y.Qlo <- quantile(frac.bounds.summary$mean - frac.bounds.summary$ci,
                         probs=c(0.1), na.rm=TRUE)[1]
    p1.y.Qhi <- quantile(frac.bounds.summary$mean + frac.bounds.summary$ci,
                         probs=c(0.9), na.rm=TRUE)[1]
    p1.y.min <- max(min(frac.bounds.summary$mean - frac.bounds.summary$ci),
                    p1.y.Qlo - 1.5*(p1.y.Qhi - p1.y.Qlo))
    p1.y.max <- min(max(frac.bounds.summary$mean + frac.bounds.summary$ci),
                    p1.y.Qhi + 1.5*(p1.y.Qhi - p1.y.Qlo))
    if (in.taxonomy %in% c("tax-11620_S", "tax-138948_None")) {
        # In main text figure; manually make clean axis breaks
        p1 <- p1 + scale_y_continuous(limits=c(p1.y.min, p1.y.max),
                                      breaks=c(40, 60, 80, 100))
    } else {
        p1 <- p1 + scale_y_continuous(limits=c(p1.y.min, p1.y.max))
    }
    # Add title to plot and axis labels
    p1 <- p1 + xlab("Genome position") + ylab("Detected sequences (%)")
    p1 <- p1 + ggtitle(taxonomy.name)
    p1 <- p1 + theme_pubr()

    # Second, produce a plot for the real.min-guides design, showing the
    # number of guides required in each window for each choice of guide
    # coverage (gp)
    # Produce plot where x-axis shows the reference position of the middle of
    # the window
    p2 <- ggplot(real.min.dist.summary, aes(x=ref.pos))
    # Plot mean values as points
    p2 <- p2 + geom_line(aes(y=mean, color=gp.display), size=1.0)
    # Can use geom_errorbar(..) to show error bars at each plotted
    # x-value; alternatively, geom_ribbon(..) to show a continuous
    # interval (i.e., confidence band) around each line. Note that
    # this is a 95% pointwise confidence band, NOT a simultaneous
    # confidence band.
    p2 <- p2 + geom_ribbon(aes(ymin=mean-ci,
                               ymax=mean+ci,
                               fill=gp.display), alpha=0.3)
    # Manually set the y-axis limits to avoid outliers in the confidence
    # intervals (the ribbon); these outliers may not show on the plot
    p2.y.Qlo <- quantile(real.min.dist.summary$mean - real.min.dist.summary$ci,
                         probs=c(0.1), na.rm=TRUE)[1]
    p2.y.Qhi <- quantile(real.min.dist.summary$mean + real.min.dist.summary$ci,
                         probs=c(0.9), na.rm=TRUE)[1]
    p2.y.min <- max(min(real.min.dist.summary$mean - real.min.dist.summary$ci),
                    p2.y.Qlo - 1.5*(p2.y.Qhi - p2.y.Qlo))
    p2.y.max <- min(max(real.min.dist.summary$mean + real.min.dist.summary$ci),
                    p2.y.Qhi + 1.5*(p2.y.Qhi - p2.y.Qlo))
    p2.y.min <- min(0, p2.y.min)    # make sure to show y=0 guides
    if (in.taxonomy == "tax-11620_S") {
        # In main text figure; manually make clean axis breaks
        p2 <- p2 + scale_y_continuous(limits=c(1, p2.y.max),
                                      breaks=c(1, 4, 7, 10))
    } else {
        p2 <- p2 + scale_y_continuous(limits=c(p2.y.min, p2.y.max))
    }
    # Add axis labels
    p2 <- p2 + xlab("Genome position") + ylab("Number of probes")
    p2 <- p2 + ggtitle(taxonomy.name)
    p2 <- p2 + theme_pubr()

    # Third, produce a plot focused on a baseline using the N most
    # common subsequences from N=1..10 (loosely called 'mode')
    # This shows the fraction of sequences detected, comparing to ADAPT
    # with 1, 2, and 3 guides
    # Because it can be hard to visualize up to N=10, this uses a stacked
    # area chart for the 'mode' baseline where each group corresponds to
    # a choice of N
    frac.bounds.summary.all.mode.adapt <- frac.bounds.summary.all.mode[frac.bounds.summary.all.mode$approach %in% paste0("real.max.activity.frac.bound.hgc-", c(1:3)),]
    frac.bounds.summary.all.mode.modes <- frac.bounds.summary.all.mode[frac.bounds.summary.all.mode$approach %in% paste0("naive.mode.frac.bound.upto-", c(1:10)),]
    # Compute how much each additional probe adds to the fraction covered
    t <- frac.bounds.summary.all.mode.modes
    t$num.probes <- as.numeric(gsub(".+upto-(\\d+)", "\\1", t$approach))
    t$delta.frac.bound <- NA
    t$delta.frac.bound[t$num.probes == 1] <- t$mean[t$num.probes == 1]  # at 1 probe, delta is against 0
    for(n in 2:10) {
        t$delta.frac.bound[t$num.probes == n] <- t$mean[t$num.probes == n] - t$mean[t$num.probes == n-1]
    }
    t$num.probes.factor <- factor(t$num.probes, levels=rev(unique(t$num.probes)))   # specify order of areas, with N=1 probe on bottom
    frac.bounds.summary.all.mode.modes <- t
    p3 <- ggplot(frac.bounds.summary.all.mode, aes(x=ref.pos))
    # Add a stacked area chart for the 'mode' baselines; show the delta in
    # fraction bound (N probes versus N-1 probes) so each area represents
    # the coverage added by adding a probe (the top area is the coverage
    # using all N=10 probes)
    # Note that default for geom_area() is stacked
    # This sets alpha according to num.probes
    p3 <- p3 + geom_area(data=frac.bounds.summary.all.mode.modes,
                         aes(x=ref.pos, y=delta.frac.bound, group=num.probes.factor, alpha=num.probes),
                         fill="black",
                         color="white", size=0.2) # thin white for outline between areas
    # Plot mean values as a line (only for ADAPT)
    p3 <- p3 + geom_line(data=frac.bounds.summary.all.mode.adapt,
                         aes(y=mean, color=approach.display), size=0.9)
    # Can use geom_errorbar(..) to show error bars at each plotted
    # x-value; alternatively, geom_ribbon(..) to show a continuous
    # interval (i.e., confidence band) around each line. Note that
    # this is a 95% pointwise confidence band, NOT a simultaneous
    # confidence band. (And note it is only for ADAPT)
    # For simplicity on this plot, leave it out
    #p3 <- p3 + geom_ribbon(data=frac.bounds.summary.all.mode.adapt,
    #                       aes(ymin=mean-ci,
    #                           ymax=mean+ci,
    #                           fill=approach.display), alpha=0.2)
    # Add title to plot and axis labels
    p3 <- p3 + xlab("Genome position") + ylab("Detected sequences (%)")
    p3 <- p3 + ggtitle(taxonomy.name)
    p3 <- p3 + theme_pubr()

    # Make sure p1 and p2 are on the same x-axis scale (same limit)
    # Force x.min to 0 even if the first window (after filtering out
    # windows at the ends) is not at 0
    # Put p3 on these limits as well
    #x.min <- min(min(frac.bounds.summary$window.start),
    #             min(real.min.dist.summary$window.start))
    x.min <- 0
    x.max <- max(max(frac.bounds.summary$window.start),
                 max(real.min.dist.summary$window.start))
    p1 <- p1 + scale_x_continuous(limits=c(x.min, x.max))
    p2 <- p2 + scale_x_continuous(limits=c(x.min, x.max))
    p3 <- p3 + scale_x_continuous(limits=c(x.min, x.max))

    # Use name="" to avoid legend title
    # Use viridis color scheme, but specified manually to draw contrast between
    # the naive designs and real designs (only important for p1)
    p1.cols <- c("Consensus"="#AC3876",
                 "Mode"="#D69856",
                 "ADAPT, 1 probe"="#3D0C51",
                 "ADAPT, 2 probes"="#5EB47E",
                 "ADAPT, 3 probes"="#FAE655")
    p2.cols <- c("90%"="#3D0C51",
                 "95%"="#5EB47E",
                 "99%"="#FAE655")
    p1 <- p1 + scale_color_manual(name="",
                                  values=p1.cols)
    p1 <- p1 + scale_fill_manual(name="",
                                 values=p1.cols)
    p2 <- p2 + scale_color_manual(name="",
                                  values=p2.cols)
    p2 <- p2 + scale_fill_manual(name="",
                                 values=p2.cols)
    p3 <- p3 + scale_color_manual(name="",
                                  values=p1.cols)
    p3 <- p3 + scale_fill_manual(name="",
                                 values=p1.cols)
    p3 <- p3 + scale_alpha_continuous(name="# subsequences",
                                      breaks=c(1,10),
                                      range=c(0.1, 0.9))

    r <- list(max.frac.bound=p1, min.num.guides=p2, max.frac.bound.modes=p3)
    return(r)
}


r <- plot.results.for.taxonomy(in.taxonomy)
height <- 3.2
width <- 2.4 * height
ggsave(file.path(out.dir, paste0(in.taxonomy, ".max-frac-bound.pdf")), r$max.frac.bound, width=width, height=height, useDingbats=FALSE)
ggsave(file.path(out.dir, paste0(in.taxonomy, ".min-num-guides.pdf")), r$min.num.guides, width=width, height=height, useDingbats=FALSE)
ggsave(file.path(out.dir, paste0(in.taxonomy, ".max-frac-bound-modes.pdf")), r$max.frac.bound.modes, width=width, height=height, useDingbats=FALSE)
