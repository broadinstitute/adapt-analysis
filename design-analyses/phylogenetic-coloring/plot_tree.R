#!/usr/bin/env Rscript

# Plot a phylogenetic tree with predicted assay performance.
#
# Args:
#  1: taxon to use (e.g., 'lasv-L')
#  2: design id (ranking) to use (e.g., '1', '2', ...)
#
# By Hayden Metsky <hmetsky@broadinstitute.org>

library(ggplot2)
library(ggtree)
require(dplyr)


args <- commandArgs(trailingOnly=TRUE)
taxon <- args[1]
assay.option <- args[2]

data.dir <- paste0("data/", taxon, "/")
out.dir <- paste0("out/", taxon, "/")

# Read previous analyses
per.seq.predictions <- data.frame(read.table(file.path(paste0(out.dir, "per-seq-predictions.tsv")), sep="\t", header=TRUE))
per.seq.metadata <- data.frame(read.table(file.path(paste0(out.dir, "per-seq-metadata.tsv")), sep="\t", header=TRUE))

# Replace '_' in column names with '.'
names(per.seq.predictions) <- gsub("_", ".", names(per.seq.predictions))
names(per.seq.metadata) <- gsub("_", ".", names(per.seq.metadata))

# Remove sequence version ([accession].[version]) from sequence names in predictions to
# get only accession
per.seq.predictions$accession <- gsub("(.+)\\.(.+)", "\\1", per.seq.predictions$seq.name)

# Pull out only the requested design assay option
per.seq.predictions <- per.seq.predictions[per.seq.predictions$design.id == assay.option,]

# Convert the ideal guide sequence for each target into an alphabetic label
# If there are two guides, the label will be 'A' or 'B' corresponding to
#   which is the best guide for each target sequence
# Note this must be done *after* pulling out the requested design assay
#   option (above), or else there will be more unique guides than we need
#   to plot (use `factor()` again to remove factors no longer used because they
#   were from other design options)
ideal.guide <- per.seq.predictions$guide.ideal.target.sequence
ideal.guide[ideal.guide == "None"] <- NA
ideal.guide.numeric <- as.numeric(factor(ideal.guide)) # '1', '2', etc. to label a guide
ideal.guide.alph <- intToUtf8(utf8ToInt('A') + ideal.guide.numeric - 1, multiple=TRUE) # 'A', 'B', etc. to label a guide
ideal.guide.alph[is.na(ideal.guide.alph)] <- "N/A"
per.seq.predictions$ideal.guide.label <- factor(ideal.guide.alph)

# Remove sub-country information (e.g., 'Nigeria: Ebonyi State') that comes after ':'
# In other words, pull out what comes before the ':'
per.seq.metadata$country <- gsub("^([^:]+)(:.*)?$", "\\1", per.seq.metadata$country)
per.seq.metadata$country[per.seq.metadata$country == "None"] <- "Unknown"
per.seq.metadata$country <- factor(per.seq.metadata$country)

# Join on accession
per.seq <- merge(x=per.seq.predictions, y=per.seq.metadata, by="accession")

# Rename seq.name to label, to work with ggtree and merge correctly, and
# move it to the first column (following https://github.com/YuLab-SMU/ggtree/issues/199#issuecomment-428528738)
# Note that naming it 'node' is bad and led to an issue downstream
#   ggtree(tree)$data has a 'node' column that is an integer, and differs
#   from the tip label. If per.seq has a 'node' column, it will join with
#   ggtree(tree) on 'node' (https://github.com/YuLab-SMU/ggtree/blob/4a86238b3000d266147bda244e5ea30193d4ad15/R/operator.R#L86)
#   and not be able to join correctly. As long as the first column corresponds
#   to the tip label (and is *not* named 'node'), regardless of what it is
#   named, will lead per.seq to be joined correctly with ggtree(tree) by
#   joining on ggtree(tree)$data's 'label' column (https://github.com/YuLab-SMU/ggtree/blob/4a86238b3000d266147bda244e5ea30193d4ad15/R/operator.R#L91)
names(per.seq)[names(per.seq) == "seq.name"] <- "label"
per.seq <- per.seq %>% relocate(label)
per.seq <- as.data.frame(per.seq)
per.seq$label <- as.character(per.seq$label)

# Read the tree file (in Newick format)
nwk <- paste0(data.dir, "aln.tree")
tree <- read.tree(nwk)

# Use `%<+%` operator in ggtree to combine tree with data
#  (defined at: https://github.com/YuLab-SMU/ggtree/blob/4a86238b3000d266147bda244e5ea30193d4ad15/R/operator.R#L51)
#  see note above about how these are joined together
p1 <- ggtree(tree) %<+% per.seq +
    geom_tippoint(aes(color=country), size=0.5) +
    geom_tiplab(aes(label=accession, color=country), align=TRUE, linesize=0.1, 
                 size=0.8, offset=0.1, hjust=0) +
    geom_tiplab(aes(label=year, color=country), align=TRUE, linetype=NA,
                size=0.8, offset=0.2, hjust=0) +
    scale_color_viridis_d(option="turbo", name="Country") +
    guides(color=guide_legend(override.aes = list(size=5))) # larger legend size

# Add heatmap of guide activity
df.guide.activity <- per.seq[c("guide.activity")]
df.guide.activity <- as.data.frame(df.guide.activity)
df.guide.activity$guide.activity <- as.numeric(as.character(df.guide.activity$guide.activity))
df.guide.activity$guide.activity[is.na(df.guide.activity$guide.activity)] <- 0
rownames(df.guide.activity) <- per.seq$label
names(df.guide.activity)[names(df.guide.activity) == "guide.activity"] <- "Guide"
p2 <- gheatmap(p1, df.guide.activity, width=.4, offset=.3, colnames=TRUE,
        colnames_offset_y=-5) + #%>% 
      #scale_x_ggtree
      scale_fill_viridis_c(option="viridis", name="Predicted\nguide activity", breaks=c(0,1,2,3), limits=c(0,3.5))

# Start a new color scale, based on: http://yulab-smu.top/treedata-book/chapter7.html#gheatmap-ggnewscale
require(ggnewscale)
p3 <- p2 + new_scale_fill()

# Add heatmap of the particular guide that does best against each
# sequence (e.g., if there are 2 guides, this heatmap will show either
# 'A' or 'B' for each guide)
df.ideal.guide <- per.seq[c("ideal.guide.label")]
df.ideal.guide <- as.data.frame(df.ideal.guide)
rownames(df.ideal.guide) <- per.seq$label
names(df.ideal.guide)[names(df.ideal.guide) == "ideal.guide.label"] <- "Ideal\nguide"
p3 <- gheatmap(p3, df.ideal.guide, width=.1, offset=1.1, colnames=TRUE,
        colnames_offset_y=-9) +
    scale_fill_viridis_d(option="mako", name="Ideal\nguide")

# Start a new color scale, based on: http://yulab-smu.top/treedata-book/chapter7.html#gheatmap-ggnewscale
require(ggnewscale)
p4 <- p3 + new_scale_fill()

# Add heatmap of primer mismatches
df.primer.mismatches <- per.seq[c("left.primer.mismatches", "right.primer.mismatches")]
df.primer.mismatches <- as.data.frame(df.primer.mismatches)
df.primer.mismatches$left.primer.mismatches[df.primer.mismatches$left.primer.mismatches == "None"] <- "0"
df.primer.mismatches$right.primer.mismatches[df.primer.mismatches$right.primer.mismatches == "None"] <- "0"
df.primer.mismatches$left.primer.mismatches <- as.factor(as.numeric(as.character(df.primer.mismatches$left.primer.mismatches)))
df.primer.mismatches$right.primer.mismatches <- as.factor(as.numeric(as.character(df.primer.mismatches$right.primer.mismatches)))
rownames(df.primer.mismatches) <- per.seq$label
names(df.primer.mismatches)[names(df.primer.mismatches) == "left.primer.mismatches"] <- "5'"
names(df.primer.mismatches)[names(df.primer.mismatches) == "right.primer.mismatches"] <- "3'"
p4 <- gheatmap(p4, df.primer.mismatches, width=.4, offset=1.8, colnames=TRUE,
        colnames_offset_y=-5) +
    scale_fill_viridis_d(option="magma", name="Primer\nmismatches")

ggsave(paste0(out.dir, "plot--assay-option-", assay.option, ".pdf"), p4, width=8, height=8, useDingbats=FALSE)
