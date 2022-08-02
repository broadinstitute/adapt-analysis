# Plot data growth over time.
#
# By Hayden Metsky <hayden@mit.edu>

require(ggplot2)
require(gridExtra)
require(dplyr)
require(viridis)
require(scales)
require(ggpubr)

IN.TABLE.WITH.SARS <- "data/counts-cumulative.with-sars.tsv.gz"
IN.TABLE.WITHOUT.SARS <- "data/counts-cumulative.no-sars.tsv.gz"
OUT.GENOMES.PDF.PREFIX.WITH.SARS <- "plots/cumulative-growth.num-genomes.with-sars"
OUT.GENOMES.PDF.PREFIX.WITHOUT.SARS <- "plots/cumulative-growth.num-genomes.no-sars"
OUT.KMERS.PDF.PREFIX <- "plots/cumulative-growth.num-kmers"


make.plots <- function(in.table, cumulative.count.colname, out.pdf.prefix, y.label) {
    # Make plots for number of genomes or number of k-mers.
    #
    # Args:
    #     in.table: path to file to read
    #     count.colname: column name in counts giving the value to compute
    #     out.pdf.prefix: path to output pdf (prefix of filename)
    #     y.label: y-axis label

    # Read table and replace '_' in column names with '.'
    # Use quote="\"" to handle apostrophe in a species name
    counts <- read.table(gzfile(in.table), header=TRUE, sep="\t",
                         quote="\"")
    names(counts) <- gsub("_", ".", names(counts))

    # Only use 2006--2022
    counts <- counts[counts$year >= 2006 & counts$year <= 2022, ]

    # Make a 'cumulative.count' column, copying cumulative.count.colname
    counts$cumulative.count <- counts[, cumulative.count.colname]

    # Make the same as above, without Influenza A/B and Rotavirus
    `%notin%` <- Negate(`%in%`)
    counts.no.flu.or.rota <- counts[counts$taxid %notin% c(11320,11520,28875), ]

    # Sum counts across all taxa for each year
    counts.sum.over.taxa <- counts %>% 
        group_by(year) %>% 
        summarize(cumulative.count=sum(cumulative.count))
    counts.sum.over.taxa.no.flu.or.rota <- counts.no.flu.or.rota %>% 
        group_by(year) %>% 
        summarize(cumulative.count=sum(cumulative.count))

    # Plot geom_area of all counts
    p <- ggplot(counts.sum.over.taxa, aes(x=year, y=cumulative.count))
    p <- p + geom_area()
    p <- p + ggtitle("Cumulative counts - all")
    p <- p + xlab("Year") + ylab(y.label)
    p <- p + scale_x_continuous(breaks=c(2006,2008,2010,2012,2014,2016,2018,2020,2022))   # use more breaks on x-axis
    p <- p + scale_y_continuous(labels=comma)   # don't use scientific notation on y-axis
    p <- p + theme_pubr()
    p.all <- p

    # Make the same as above, without Influenza A/B and Rotavirus
    p <- ggplot(counts.sum.over.taxa.no.flu.or.rota, aes(x=year, y=cumulative.count))
    p <- p + geom_area()
    p <- p + ggtitle("Cumulative counts - without Influenza A/B or Rotavirus A")
    p <- p + xlab("Year") + ylab(y.label)
    p <- p + scale_x_continuous(breaks=c(2006,2008,2010,2012,2014,2016,2018,2020,2022))   # use more breaks on x-axis
    p <- p + scale_y_continuous(labels=comma)   # don't use scientific notation on y-axis
    p <- p + theme_pubr()
    p.no.flu.or.rota <- p

    # Plot geom_area of all counts, stacked by taxa
    counts.last.year <- counts[counts$year == 2018, ]
    taxa.order <- counts.last.year[order(counts.last.year$cumulative.count), ]$taxid
    counts.cp <- data.frame(counts)
    counts.cp$taxid <- factor(as.character(counts.cp$taxid), levels=taxa.order)
    p <- ggplot(counts.cp, aes(x=year, y=cumulative.count))
    p <- p + geom_area(aes(fill=taxid))
    p <- p + ggtitle("Cumulative counts - all")
    p <- p + xlab("Year") + ylab(y.label)
    p <- p + scale_x_continuous(breaks=c(2006,2008,2010,2012,2014,2016,2018,2020,2022))   # use more breaks on x-axis
    p <- p + scale_y_continuous(labels=comma)   # don't use scientific notation on y-axis
    # If we use a regular scale, it will use a number of colors equal to the
    # number of taxid and gradually change between these, making it impossible
    # to tell apart adjacent species in the stacked area
    # Instead pick out 5 colors and alternate between these
    colors <- magma(5)  # color scale from viridis package
    p <- p + scale_fill_manual(values=rep_len(colors, length(unique(counts$taxid))))
    p <- p + guides(fill=FALSE) # no legend
    p <- p + theme_pubr()
    p.stacked <- p

    # Save PDFs
    ggsave(paste0(out.pdf.prefix, ".all.pdf"), p.all, width=5, height=5, useDingbats=FALSE)
    ggsave(paste0(out.pdf.prefix, ".no-flu-or-rota.pdf"), p.no.flu.or.rota, width=5, height=5, useDingbats=FALSE)
    ggsave(paste0(out.pdf.prefix, ".stacked.pdf"), p.stacked, width=5, height=5, useDingbats=FALSE)

    # Delete the column created at the start
    counts$cumulative.count <- NULL
}

make.plots(IN.TABLE.WITH.SARS, "cumulative.num.genomes", OUT.GENOMES.PDF.PREFIX.WITH.SARS, "Cumulative number of genomes")
make.plots(IN.TABLE.WITHOUT.SARS, "cumulative.num.genomes", OUT.GENOMES.PDF.PREFIX.WITHOUT.SARS, "Cumulative number of genomes")
#make.plots(IN.TABLE.WITH.SARS, "cumulative.num.unique.kmers", OUT.KMERS.PDF.PREFIX, "Cumulative number of unique k-mers")
