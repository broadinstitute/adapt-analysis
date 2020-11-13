# Plot variants along genome at different time points.
#
# This plots the genome on the x-axis and time on the y-axis, going
# downward. At each date, it plots variants *added* as called in
# genomes from that date compared to the previous one; the only variants
# shown are ones also called in the genomes at the final date.
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
last.date <- levels(variants$date)[[length(levels(variants$date))]]

max.genome.pos <- max(variants$pos)

# Do not keep variants in the first or last 100 nt; there may be many variants
# here in part due to assembly, alignment, or sequencing artifacts
# The first and last 100 nt are still within the 5' and 3' UTRs, respectively,
# for SARS-CoV-2
variants <- variants[variants$pos >= 100 & variants$pos <= max.genome.pos - 100,]

# Call variant sites as ones where the fraction of genomes counted in which there is
# a variant is >= 0.1%; this makes sense when the input is *not* subsampled
# If the input were subsampled to a fixed number for each date, then we may
# want to call variants based on when the number exceeds a fixed number (like
# 2, to avoid singletons)
variants$freq <- variants$num.genomes.with.variant / variants$num.genomes.counted
variants <- variants[variants$freq >= 0.001,]
#variants <- variants[variants$num.genomes.with.variant >= 2,]

# Set whether each variant is low (<1%) or high (>=1%) frequency, to adjust
# its color
variants$freq.level <- ifelse(variants$freq < 0.01, "low", "high")

# Label every variant as either being present in the first date (at high
# frequency) as ancestral, or arising since then (new)
variants$status <- "new"
ancestral.high.variant.pos <- variants[variants$date == first.date & variants$freq.level == "high",]$pos
variants[variants$pos %in% ancestral.high.variant.pos,]$status <- "ancestral.high"
ancestral.low.variant.pos <- variants[variants$date == first.date & variants$freq.level == "low",]$pos
variants[variants$pos %in% ancestral.low.variant.pos,]$status <- "ancestral.low"

# Make a dataframe that has added variants at every date -- i.e., ones that
# will be shown in last.date, and it only keeps the row when that variant
# first got called (or 'captured' in the data)
# Only show variants at the frequency level they hit on last.date -- i.e., if
# a variant is high frequency in last.date, remove it from added.variants
# when it is low frequency, so it is only shown as added when it first is
# at high frequency
# For variants to show on the last date: only keep ones that either: (a) are
# not at all ancestral; or (b) are low frequency in ancestral and then rose to
# high frequency (do *not* show ones that are low frequency in ancestral and
# stayed low frequency or are high frequency in ancestral)
variants.only.last.date <- variants[variants$date == last.date & variants$status != "ancestral.high" & (variants$status != "ancestral.low" | variants$freq.level == "high"),]
variants.only.last.date$final.freq <- variants.only.last.date$freq
variants.only.last.date$final.freq.level <- ifelse(variants.only.last.date$freq < 0.01, "low", "high")
added.variants <- variants[variants$pos %in% variants.only.last.date$pos & variants$date != first.date,]
added.variants$final.freq <- with(variants.only.last.date, freq[match(added.variants$pos, pos)])
added.variants$final.freq.level <- ifelse(added.variants$final.freq < 0.01, "low", "high")
added.variants <- added.variants[added.variants$freq.level == added.variants$final.freq.level,]
# De-duplicating will select the first of ambiguous columns, which is this case
# will correspond to the first date for each added variant
added.variants <- added.variants[!duplicated(added.variants$pos),]

# Print number of variants at the end
print(paste("Total number of variants in Combined:", nrow(variants.only.last.date)))
print(paste("Number of low frequency variants in Combined:", sum(variants.only.last.date$final.freq.level == "low")))
print(paste("Number of high frequency variants in Combined:", sum(variants.only.last.date$final.freq.level == "high")))

# Make a data frame to plot that combines added.variants with
# variants.only.last.date
variants.only.last.date$is.combined <- TRUE
added.variants$is.combined <- FALSE
plot.df <- rbind(added.variants, variants.only.last.date)

# Reverse the levels in date so that the most recent date is on bottom rather
# than on top
# We could do this directly in the plot with
# `scale_y_discrete(limits = rev(levels(df$date.simple)))`, but this
# will not work well with geom_segment for the variants
plot.df$date <- factor(plot.df$date, levels=rev(levels(plot.df$date)))

# Add column giving month/day only -- e.g., 'Feb. 1' instead of '2020-02-01'
# Keep the levels sorted in the right order (i.e., as given in date)
plot.df$date.simple <- factor(plot.df$date)
levels(plot.df$date.simple) <- format(as.Date(levels(plot.df$date.simple)), "%b %e")

# Add a y.label, which is date.simple for all the added.variants and "Combined"
# for the combined variants
plot.df$y.label <- factor(plot.df$date.simple, levels=c("Combined", levels(plot.df$date.simple)))
plot.df[plot.df$is.combined == TRUE,]$y.label <- "Combined"

# Have horizontal bars to represent genomes; make separate data frame for
# these genomes, with one row for each
genome.bars <- data.frame(y=levels(plot.df$y.label))
genome.bars$y <- factor(genome.bars$y, levels=levels(plot.df$y.label))
genome.bars$x.start <- rep(0, nrow(genome.bars))
genome.bars$x.end <- rep(max.genome.pos, nrow(genome.bars))

# Since data is plotted in the same order as the data frame, put the
# low frequency variants first so the high frequency ones are plotted on top
# if they overlap
plot.df$freq.level <- factor(plot.df$freq.level, levels=c("low", "high"))
plot.df <- plot.df[order(plot.df$freq.level),]

# Make a plot
# There are few ways to plot variants:
#  - geom_segment to put rectangles at variant sites
#  - geom_point, although there is no rectangular point shape
#  - geom_tile, which could work but may be complicated if we want separation
#    between genomes on the vertical axis
p <- ggplot(plot.df, aes(x=pos, y=y.label)) +
        geom_rect(data=genome.bars, inherit.aes=FALSE,   # genome backgrounds
                   aes(xmin=x.start, xmax=x.end,
                       ymin=as.numeric(y)-0.4, ymax=as.numeric(y)+0.4),
                   alpha=0.05, fill="black") +
        geom_segment(aes(xend=pos,  # variants
                         y=as.numeric(y.label)-0.4, yend=as.numeric(y.label)+0.4,
                         color=freq.level),
                     size=0.7) +
        scale_y_discrete(limits = levels(plot.df$y.label)) +  # show dates on y-axis
        scale_color_manual(guide=FALSE, # no legend
                           values=c("low"="#CCBAD1", "high"="#430053")) +   # two colors from viridis scale
        theme_pubr() +
        xlab("Genome") +
        ylab("Date")

ggsave(OUT.PDF, p, device="pdf", width=6, height=3.5, useDingbats=FALSE)
