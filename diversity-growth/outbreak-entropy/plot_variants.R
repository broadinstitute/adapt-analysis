# Plot variants along genome at different time points.
#
# This plots the genome on the x-axis and time on the y-axis, going
# downward.
#
# By Hayden Metsky <hayden@mit.edu>


require(ggplot2)
require(viridis)
require(ggpubr)

args <- commandArgs(trailingOnly=TRUE)
IN.TSV <- args[1]
OUT.PDF <- args[2]


# Read TSV of results and replace '_' in column names with '.'
variants <- read.table(gzfile(IN.TSV), header=TRUE, sep="\t")
names(variants) <- gsub("_", ".", names(variants))

first.date <- levels(variants$date)[[1]]

# Reverse the levels in date so that the most recent date is on bottom rather
# than on top
# We could do this directly in the plot with
# `scale_y_discrete(limits = rev(levels(variants$date.simple)))`, but this
# will not work well with geom_segment for the variants
variants$date <- factor(variants$date, levels=rev(levels(variants$date)))

# Add column giving month/day only -- e.g., 'Feb. 1' instead of '2020-02-01'
# Keep the levels sorted in the right order (i.e., as given in date)
variants$date.simple <- factor(variants$date)
levels(variants$date.simple) <- format(as.Date(levels(variants$date.simple)), "%b %e")

max.genome.pos <- max(variants$pos)

# Do not keep variants in the first or last 100 nt; there may be many variants
# here in part due to assembly, alignment, or sequencing artifacts
# The first and last 100 nt are still within the 5' and 3' UTRs, respectively,
# for SARS-CoV-2
variants <- variants[variants$pos >= 100 & variants$pos <= max.genome.pos - 100,]

# Call variant sites as ones where the fraction of genomes counted in which there is
# a variant is >= 1%; this makes sense when the input is *not* subsampled
# If the input were subsampled to a fixed number for each date, then we may
# want to call variants based on when the number exceeds a fixed number (like
# 2, to avoid singletons)
variants <- variants[variants$num.genomes.with.variant / variants$num.genomes.counted >= 0.01,]
#variants <- variants[variants$num.genomes.with.variant >= 2,]

# Label every variant as either being present in the first date, or arising
# since then ('ancestral' or 'new')
variants$status <- "new"
ancestral.variant.pos <- variants[variants$date == first.date,]$pos
variants[variants$pos %in% ancestral.variant.pos,]$status <- "ancestral"

# Have horizontal bars to represent genomes; make separate data frame for
# these genomes, with one row for each
genome.bars <- data.frame(y=levels(variants$date.simple))
genome.bars$y <- factor(genome.bars$y, levels=levels(variants$date.simple))
genome.bars$x.start <- rep(0, nrow(genome.bars))
genome.bars$x.end <- rep(max.genome.pos, nrow(genome.bars))

# Make a plot
# There are few ways to plot variants:
#  - geom_segment to put rectangles at variant sites
#  - geom_point, although there is no rectangular point shape
#  - geom_tile, which could work but may be complicated if we want separation
#    between genomes on the vertical axis
p <- ggplot(variants, aes(x=pos, y=date.simple)) +
        geom_rect(data=genome.bars, inherit.aes=FALSE,   # genome backgrounds
                   aes(xmin=x.start, xmax=x.end,
                       ymin=as.numeric(y)-0.35, ymax=as.numeric(y)+0.35),
                   alpha=0.07, fill="black") +
        geom_segment(aes(xend=pos,
                         y=as.numeric(date.simple)-0.35, yend=as.numeric(date.simple)+0.35,
                         color=status),
                     size=1.5) +
        scale_y_discrete(limits = levels(variants$date.simple)) +  # show dates on y-axis
        scale_color_manual(guide=FALSE, # no legend
                           values=c("ancestral"="#60CA68", "new"="#430053")) +   # two colors from viridis scale
        theme_pubr() +
        xlab("Genome") +
        ylab("Date")

ggsave(OUT.PDF, p, device="pdf", width=16, height=8, useDingbats=FALSE)
