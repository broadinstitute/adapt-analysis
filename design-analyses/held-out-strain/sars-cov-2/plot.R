#!/usr/bin/env Rscript

# Plot results of analyzing coverage on SARS-CoV-2 for
# designs made in 2018.
#
# By Hayden Metsky <hmetsky@broadinstitute.org>

library(ggplot2)
library(gridExtra)
library(reshape2)
library(viridis)
library(ggpubr)

IN.FN <- "evaluations/evaluations-compiled.tsv"
vals <- data.frame(read.table(IN.FN, header=TRUE, quote="", sep="\t"))

# Replace '_' in column names with '.'
names(vals) <- gsub("_", ".", names(vals))

plot.stat <- function(stat, design, y.label, out.pdf) {
    # Args:
    #   stat: name of statistic in input file
    #   design: name of design in input file
    #   y.label: y-axis label
    #   out.pdf: path to output PDF

    # Pull out the statistic to plot
    vals.stat <- vals[vals$stat == stat,]

    # Pull out the design to plot
    vals.stat <- vals.stat[vals.stat$design == design,]

    # Only keep the top 5 designs
    vals.stat <- vals.stat[vals.stat$design.id <= 5,]

    # Make design.id a factor
    vals.stat$design.id <- factor(vals.stat$design.id)

    # Make plot of value for each design, separated by whether
    # it was evaluated against all input or SARS-CoV-2
    p <- ggplot(vals.stat, aes(x=design.id, y=value)) +
        geom_point(aes(color=evaluated.against),
                   size=3, position=position_dodge(width=0.3)) +
        xlab("Design ranking") + ylab(y.label) +
        scale_color_manual(values=c("input"="#643C72", "sars-cov-2"="#7DBE9C"),
                           labels=c("input"="Design input", "sars-cov-2"="SARS-CoV-2")) +
        theme_pubr() +
        theme(legend.title=element_blank())

    # Save PDF
    ggsave(out.pdf, p, width=4, height=4, useDingbats=FALSE)
}

plot.stat("frac-bound-active", "all", "Fraction of genomes detected", "plots/all.frac-bound-active.pdf")
plot.stat("mean-guide-activity", "all", "Mean activity of guide set", "plots/all.mean-guide-activity.pdf")
plot.stat("frac-bound-active", "downsampled-sars-cov-1", "Fraction of genomes detected", "plots/downsampled-sars-cov-1.frac-bound-active.pdf")
plot.stat("mean-guide-activity", "downsampled-sars-cov-1", "Mean activity of guide set", "plots/downsampled-sars-cov-1.mean-guide-activity.pdf")
