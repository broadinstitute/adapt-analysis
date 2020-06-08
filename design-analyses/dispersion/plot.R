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
library(ggpubr)


# Manually specify taxonomies
TAXONOMIES <- c("tax-64320_None",
                "tax-121791_None",
                "tax-186538_None",
                "tax-11320_2",
                "tax-11620_S",
                "tax-138948_None")
TAXONOMY.NAMES <- c("Zika",
                    "Nipah",
                    "Ebola",
                    "IAV",
                    "Lassa",
                    "Enterovirus A")


args <- commandArgs(trailingOnly=TRUE)
out.dir <- args[1]


plot.violin <- function(filename, title) {
    # Args:
    #   filename: name of file containing distribution (one
    #       value per line), inside a taxonomy directory
    #   title: title of plot

    # Fill in a data frame with a distribution for each taxonomy
    dist <- data.frame(taxonomy=character(), taxonomy.name=character(), value=double())
    for (i in 1:length(TAXONOMIES)) {
        taxonomy <- TAXONOMIES[[i]]
        taxonomy.name <- TAXONOMY.NAMES[[i]]
        dist.tax <- data.frame(read.table(file.path(taxonomy, filename)))
        colnames(dist.tax) <- c("value")
        dist.tax$taxonomy <- rep(taxonomy, nrow(dist.tax))
        dist.tax$taxonomy.name <- rep(taxonomy.name, nrow(dist.tax))
        dist <- rbind(dist, dist.tax)
    }

    p <- ggplot(dist, aes(x=taxonomy.name, y=value))

    # Show a violin plot
    # Use scale="width" (instead of the default, scale="area") so
    # that all violins have the same width
    p <- p + geom_violin(trim=TRUE,
                         scale="width",
                         adjust=2,
                         fill="black",
                         alpha=0.5)

    # Add a bar with the mean +/- 1 std deviation (1 because mult=1)
    p <- p + stat_summary(fun.data="mean_sdl", fun.args=list(mult=1),
                          geom="pointrange", color="black")

    p <- p + xlab("Species") + ylab("Pairwise Jaccard similarity")

    p <- p + theme_pubr()
    p <- p + theme(axis.text.x=element_text(angle=45, hjust=1), # axis text on 45 degree angle
                   strip.background=element_blank())

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

ggsave(file.path(out.dir, "dispersion.non-resampled.exact.pdf"), p1, width=4, height=4, useDingbats=FALSE)
ggsave(file.path(out.dir, "dispersion.non-resampled.loose-equality.pdf"), p2, width=4, height=4, useDingbats=FALSE)
ggsave(file.path(out.dir, "dispersion.resampled.exact.pdf"), p3, width=4, height=4, useDingbats=FALSE)
ggsave(file.path(out.dir, "dispersion.resampled.loose-equality.pdf"), p4, width=4, height=4, useDingbats=FALSE)
