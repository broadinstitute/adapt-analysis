#!/usr/bin/env Rscript

# Plot results of cross-validation.
#
# This uses the results of analyses that were already performed.
#
# Args:
#  1: path to output PDF
#
# By Hayden Metsky <hayden@mit.edu>

library(ggplot2)
library(gridExtra)


taxonomies <- list.files(path=".", pattern="tax-*")

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

    # Multiply fractions by 100 to obtain percents
    dist$value <- dist$value * 100

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
    p <- p + xlab("Taxonomy") + ylab("Coverage (%)")

    # Leave out usual ggplot2 background and grid lines, but keep border
    # Use aspect.ratio=1 to make the plot square
    p <- p + theme_bw()
    p <- p + theme(strip.background=element_blank(),
                   panel.grid.major=element_blank(),
                   panel.grid.minor=element_blank(),
                   panel.border=element_rect(colour="black"),
                   aspect.ratio=1,
                   plot.title=element_text(size=12))

    return(p)
}


p1 <- plot.violin("coverage-against-test.distribution.txt",
                  "Cross-validation")

g <- grid.arrange(p1)
ggsave(out.pdf, g, width=8, height=8, useDingbats=FALSE)
