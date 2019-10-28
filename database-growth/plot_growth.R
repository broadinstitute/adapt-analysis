# Plot data growth over time.
#
# By Hayden Metsky <hayden@mit.edu>

require(ggplot2)
require(gridExtra)
require(dplyr)
require(viridis)
require(scales)

IN.TABLE <- "data/counts-cumulative.tsv.gz"
OUT.PDF <- "plots/cumulative-growth.pdf"


# Read table and replace '_' in column names with '.'
# Use quote="\"" to handle apostrophe in a species name
counts <- read.table(gzfile(IN.TABLE), header=TRUE, sep="\t",
                     quote="\"")
names(counts) <- gsub("_", ".", names(counts))

# Only use 2005--2018
counts <- counts[counts$year >= 2005 & counts$year <= 2018, ]

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
p <- p + xlab("Year") + ylab("Total number of genomes")
p <- p + scale_x_continuous(breaks=c(2006,2008,2010,2012,2014,2016,2018))   # use more breaks on x-axis
p <- p + scale_y_continuous(labels=comma)   # don't use scientific notation on y-axis
p <- p + theme_bw()
p.all <- p

# Make the same as above, without Influenza A/B and Rotavirus
p <- ggplot(counts.sum.over.taxa.no.flu.or.rota, aes(x=year, y=cumulative.count))
p <- p + geom_area()
p <- p + ggtitle("Cumulative counts - without Influenza A/B or Rotavirus A")
p <- p + xlab("Year") + ylab("Total number of genomes")
p <- p + scale_x_continuous(breaks=c(2006,2008,2010,2012,2014,2016,2018))   # use more breaks on x-axis
p <- p + scale_y_continuous(labels=comma)   # don't use scientific notation on y-axis
p <- p + theme_bw()
p.no.flu.or.rota <- p

# Plot geom_area of all counts, stacked by taxa
counts.last.year <- counts[counts$year == 2018, ]
taxa.order <- counts.last.year[order(counts.last.year$cumulative.count), ]$taxid
counts$taxid <- factor(as.character(counts$taxid), levels=taxa.order)
p <- ggplot(counts, aes(x=year, y=cumulative.count))
p <- p + geom_area(aes(fill=taxid))
p <- p + ggtitle("Cumulative counts - all")
p <- p + xlab("Year") + ylab("Total number of genomes")
p <- p + scale_x_continuous(breaks=c(2006,2008,2010,2012,2014,2016,2018))   # use more breaks on x-axis
p <- p + scale_y_continuous(labels=comma)   # don't use scientific notation on y-axis
# If we use a regular scale, it will use a number of colors equal to the
# number of taxid and gradually change between these, making it impossible
# to tell apart adjacent species in the stacked area
# Instead pick out 5 colors and alternate between these
colors <- magma(5)  # color scale from viridis package
p <- p + scale_fill_manual(values=rep_len(colors, length(unique(counts$taxid))))
p <- p + guides(fill=FALSE) # no legend
p <- p + theme_bw()
p.stacked <- p

# Save PDF
ggsave(OUT.PDF, arrangeGrob(p.all,
                            p.no.flu.or.rota,
                            p.stacked,
                            ncol=1),
       width=8, height=24, useDingbats=FALSE)

