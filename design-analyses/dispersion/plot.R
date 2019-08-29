#!/usr/bin/env Rscript

# Plot statistics on dispersion.
#
# This uses the results of analyses that were already performed.
#
# Args:
#  1: path to output PDF
#
# By Hayden Metsky <hayden@mit.edu>

library(ggplot2)
library(gridExtra)


taxonomies <- list.files(path=".", pattern="^tax-*")

args <- commandArgs(trailingOnly=TRUE)
out.pdf <- args[1]


plot.violin <- function(filename, title) {
    # Args:
    #   filename: name of file containing distribution (one
    #       value per line), inside a taxonomy directory
    #   title: title of plot

    # Fill in a data frame with a distribution for each taxonomy
    dist <- data.frame(taxonomy=character(), value=double())
    for (taxonomy in taxonomies) {
        dist.tax <- data.frame(read.table(file.path(taxonomy, filename)))
        colnames(dist.tax) <- c("value")
        dist.tax$taxonomy <- rep(taxonomy, nrow(dist.tax))
        dist <- rbind(dist, dist.tax)
    }

    p <- ggplot(dist, aes(x=taxonomy, y=value))

    # Show a violin plot
    # Use scale="width" (instead of the default, scale="area") so
    # that all violins have the same width
    p <- p + geom_violin(trim=TRUE,
                         scale="width",
                         adjust=2)

    # Add a bar with the mean +/- 1 std deviation (1 because mult=1)
    p <- p + stat_summary(fun.data="mean_sdl", fun.args=list(mult=1),
                          geom="pointrange", color="lightgray")

    # Add title to plot and axis labels
    p <- p + ggtitle(title)
    p <- p + xlab("Taxonomy") + ylab("Pairwise Jaccard similarity")

    # Leave out usual ggplot2 background and grid lines, but keep border
    # Use aspect.ratio=1 to make the plot square
    p <- p + theme_bw()
    p <- p + theme(axis.text.x=element_text(angle=90, hjust=1),
                   strip.background=element_blank(),
                   panel.grid.major=element_blank(),
                   panel.grid.minor=element_blank(),
                   panel.border=element_rect(colour="black"),
                   aspect.ratio=1,
                   plot.title=element_text(size=12))

    return(p)
}


p1 <- plot.violin("dispersion.non-resampled.distribution.txt",
                  "Dispersion: non-resampled")
p2 <- plot.violin("dispersion.non-resampled.loose-equality.distribution.txt",
                  "Dispersion: non-resampled (loose equality)")
p3 <- plot.violin("dispersion.resampled.distribution.txt",
                  "Dispersion: resampled")
p4 <- plot.violin("dispersion.resampled.loose-equality.distribution.txt",
                  "Dispersion: resampled (loose equality)")

g <- grid.arrange(p1, p2, p3, p4)
ggsave(out.pdf, g, width=8, height=8, useDingbats=FALSE)
