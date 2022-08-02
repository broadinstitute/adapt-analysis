# Plot data growth over time.
#
# By Hayden Metsky <hayden@mit.edu>

require(ggplot2)
require(gridExtra)
require(dplyr)
require(viridis)
require(scales)
require(ggpubr)


args <- commandArgs(trailingOnly=TRUE)
IN.TSV <- args[1]
OUT.PDF <- args[2]


make.plots <- function(in.table, cumulative.count.colname, out.pdf, y.label) {
    # Make plots for number of genomes or number of k-mers.
    #
    # Args:
    #     in.table: path to file to read
    #     count.colname: column name in counts giving the value to compute
    #     out.pdf: path to output pdf 
    #     y.label: y-axis label

    # Read table and replace '_' in column names with '.'
    # Use quote="\"" to handle apostrophe in a species name
    counts <- read.table(gzfile(in.table), header=TRUE, sep="\t",
                         quote="\"")
    names(counts) <- gsub("_", ".", names(counts))

    # Only use 2002--2022
    counts <- counts[counts$year >= 2002 & counts$year <= 2022, ]

    # Make a 'cumulative.count' column, copying cumulative.count.colname
    counts$cumulative.count <- counts[, cumulative.count.colname]

    # Plot geom_area of all counts
    p <- ggplot(counts, aes(x=year, y=cumulative.count))
    p <- p + geom_area()
    p <- p + ggtitle("Cumulative counts - all")
    p <- p + xlab("Year") + ylab(y.label)
    p <- p + scale_x_continuous(breaks=c(2002,2004,2006,2008,2010,2012,2014,2016,2018,2020,2022))   # use more breaks on x-axis
    p <- p + scale_y_continuous(labels=comma)   # don't use scientific notation on y-axis
    p <- p + theme_pubr()
    p.all <- p

    # Save PDFs
    ggsave(out.pdf, p.all, width=5, height=5, useDingbats=FALSE)

    # Delete the column created at the start
    counts$cumulative.count <- NULL
}

make.plots(IN.TSV, "cumulative.num.sequences", OUT.PDF, "Cumulative number of sequences")
