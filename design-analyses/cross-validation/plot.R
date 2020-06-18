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
library(reshape2)
library(ggpubr)
library(viridis)
library(plyr)


taxonomies <- list.files(path=".", pattern="^tax-*")

args <- commandArgs(trailingOnly=TRUE)
out.pdf <- args[1]

# Make vectors for mapping tax to names
taxs <- c("tax-11620_L",
          "tax-11620_S",
          "tax-11676_None",
          "tax-121791_None",
          "tax-138948_None",
          "tax-147711_None",
          "tax-186538_None",
          "tax-64320_None")
tax.names <- c("LASV L",
               "LASV S",
               "HIV-1",
               "NIPV",
               "EVA",
               "RVA",
               "EBOV",
               "ZIKV")


plot.violin <- function(filename, title) {
    # Args:
    #   filename: name of file containing distribution (one
    #       value per line), inside a taxonomy directory
    #   title: title of plot

    # Fill in a data frame with a distribution for each taxonomy
    dist <- data.frame(taxonomy=character(), active=double(), highly.active=double())
    for (taxonomy in taxonomies) {
        dist.tax <- data.frame(read.table(file.path(taxonomy, filename), header=TRUE))
        dist.tax$taxonomy <- rep(taxonomy, nrow(dist.tax))
        dist <- rbind(dist, dist.tax)
    }

    # Replace '_' in column names with '.'
    names(dist) <- gsub("_", ".", names(dist))

    # Sort by mean active fraction for each taxonomy (in reverse), so that
    # taxa with the highest active fraction detected come first
    dist.summarized <- ddply(dist, .(taxonomy), summarize, active=mean(active))
    dist.summarized <- dist.summarized[order(-dist.summarized$active),] # '-' for descending
    tax.ordered <- dist.summarized$taxonomy
    dist$taxonomy <- factor(dist$taxonomy, levels=tax.ordered)

    # Melt the data so active and highly.active fractions are separated
    # into different rows, rather than two columns
    dist <- melt(dist, id.vars=c("taxonomy"),
                 variable.name="threshold",
                 value.name="frac.detected")

    # Give taxa names
    dist$taxonomy.name <- mapvalues(dist$taxonomy, taxs, tax.names)

    # Multiply fractions by 100 to obtain percents
    dist$frac.detected <- dist$frac.detected * 100

    p <- ggplot(dist, aes(x=taxonomy.name, y=frac.detected, group=interaction(taxonomy.name,threshold)))

    # Specify an amount by which to space the active and highly.active
    # distributions
    dodge_width <- 0.7

    # TODO try pointrange instead of violin

    # Show a violin plot
    # Use scale="width" (instead of the default, scale="area") so
    # that all violins have the same width
    p <- p + geom_violin(aes(fill=threshold),
                         #trim=TRUE,
                         scale="width",
                         width=0.5, # shrink
                         adjust=1,
                         color=NA,  # no outline
                         position=position_dodge(width=dodge_width))

    # Add a dot for the mean
    p <- p + stat_summary(fun="mean",
                          geom="point",
                          shape=21,
                          size=2,
                          fill="black",
                          show.legend=FALSE,    # do not include in legend
                          position=position_dodge(width=dodge_width))

    # Add axis labels
    p <- p + xlab("Species") + ylab("Detected sequences (%)")

    # Leave out usual ggplot2 background and grid lines, but keep border
    # Use aspect.ratio=1 to make the plot square
    p <- p + theme_pubr()
    p <- p + theme(axis.text.x=element_text(angle=45, hjust=1)) # axis text on 45 degree angle

    # Use name="" to avoid legend title
    # Use viridis color scheme, but specified manually to avoid the yellow
    cols <- c("active"="#643C72",
              "highly.active"="#7DBE9C")
    labels <- c("active"="Active",
                "highly.active"="Highly active")
    p <- p + scale_fill_manual(name="",
                               values=cols,
                               labels=labels)

    return(p)
}


p1 <- plot.violin("coverage-against-test.distribution.txt",
                  "Cross-validation")

g <- grid.arrange(p1)
ggsave(out.pdf, g, width=8, height=8, useDingbats=FALSE)
