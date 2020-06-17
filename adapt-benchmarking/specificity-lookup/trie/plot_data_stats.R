# Plot statistics about k-mer occurrences and taxonomies.
#
# In all of these plots, k=28.
#
# Note that the k-mer occurrence distributions are only looking at *exact*
# matches across taxonomies; this may not be so useful for 28-kmers, as
# there are likely to be many more inexact matches.
#
# By Hayden Metsky <hayden@mit.edu>


require(ggplot2)
require(reshape2)
require(gridExtra)

KMER.OCC.IN.TABLE <- "out/stats.kmer-occ.tsv"
KMER.OCC.OUT.PDF <- "out/stats.kmer-occ.pdf"
TAX.STATS.IN.TABLE <- "out/stats.tax.tsv"
TAX.STATS.OUT.PDF <- "out/stats.tax.pdf"

# Read tables and replace '_' in column names with '.'
kmer.occ <- read.table(KMER.OCC.IN.TABLE, header=TRUE, sep="\t")
names(kmer.occ) <- gsub("_", ".", names(kmer.occ))

tax.stats <- read.table(TAX.STATS.IN.TABLE, header=TRUE, sep="\t")
names(tax.stats) <- gsub("_", ".", names(tax.stats))


# Compute the cumulative fraction of all k-mers that are present in
# <= X taxonomies for each X
# (Note that kmer.occ has counts for unique, umambiguous k-mers)
total.num.kmers <- sum(kmer.occ$num.kmers.in.X.taxonomies)
kmer.occ$frac.kmers.in.up.to.X.taxonomies <- cumsum(kmer.occ$num.kmers.in.X.taxonomies) / total.num.kmers


# Plot a histogram and empirical CDF of k-mer occurrences in X taxonomies
p1 <- ggplot(kmer.occ, aes(x=X, y=num.kmers.in.X.taxonomies))
p1 <- p1 + geom_bar(stat="identity")
p1 <- p1 + xlim(0, 40)
p1 <- p1 + scale_y_log10(limits=c(1,1e8))
p1 <- p1 + geom_text(aes(label=num.kmers.in.X.taxonomies), vjust=-0.25, size=1)
p1 <- p1 + xlab("X") + ylab("Number of unique k-mers present in X taxonomies")

p2 <- ggplot(kmer.occ, aes(x=X, y=frac.kmers.in.up.to.X.taxonomies))
p2 <- p2 + geom_line()
p2 <- p2 + xlim(0, 40)
p2 <- p2 + xlab("X") + ylab("Fraction of unique k-mers present in <= X taxonomies")

ggsave(KMER.OCC.OUT.PDF, arrangeGrob(p1, p2), width=8, height=8, useDingbats=FALSE)


# Plot distribution of number of k-mers in each taxonomy
p1 <- ggplot(tax.stats, aes(num.kmers.unambig))
p1 <- p1 + geom_density()
p1 <- p1 + scale_x_log10(limits=c(100,1e9))
p1 <- p1 + xlab("Number of (unambiguous) k-mers in a taxon") + ylab("Density")

p2 <- ggplot(tax.stats, aes(num.unique.kmers.unambig))
p2 <- p2 + geom_density()
p2 <- p2 + scale_x_log10(limits=c(100,1e9))
p2 <- p2 + xlab("Number of (unambiguous) unique k-mers in a taxon") + ylab("Density")

ggsave(TAX.STATS.OUT.PDF, arrangeGrob(p1, p2), width=8, height=8, useDingbats=FALSE)


# Remove the empty Rplots.pdf created above
file.remove("Rplots.pdf")
