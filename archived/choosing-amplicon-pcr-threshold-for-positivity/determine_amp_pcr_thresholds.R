#!/usr/bin/env Rscript

# Determine thresholds on amplicon pool concentration for deciding
# which samples are "positive".
#
# Args:
#  1: CSV file giving each replicate by amplicon sequencing and
#     its pool concentrations and amount sequenced
#  2: output PDF file
#  3: output TSV file for accuracy data values
#  4: output TSV file for expected num genomes data values
#
# By Hayden Metsky <hayden@mit.edu>

library(ggplot2)
library(gridExtra)

args <- commandArgs(trailingOnly=TRUE)
in.csv <- args[1]
out.pdf <- args[2]
out.accuracies.tsv <- args[3]
out.expected.nnm.gnm.tsv <- args[4]

# Stop Rplots.pdf from being made
pdf(NULL)

################
# Read the input
################
replicates <- read.table(in.csv, header=TRUE, sep=",")

# Some rows have NA for an amplicon primer pool concentration; remove these
replicates <- na.omit(replicates)

replicates$min.pool <- pmin(replicates$tapestation.primer.pool.1,
                            replicates$tapestation.primer.pool.2)
replicates$frac.covered <- replicates$unambig.bases / 10676.0

############################
# Define objective functions
############################
total.frac.covered <- sum(replicates$frac.covered)

accuracy.measure <- function(min.pool.thres) {
    # Calculate F1 score and other measures of accuracy given
    # a threshold, min.pool.thres, on the amplicon pool
    # concentration
    # Any replicates with min.pool >= min.pool.thres are labeled
    # as "positive"
    # Returns (min.pool.thres, precision, recall, F_1 score, F_0.5 score)

    # Pull out the replicates that are labeled as positive
    replicates.labeled.positive <- replicates[replicates$min.pool >= min.pool.thres, ]
    num.positive <- nrow(replicates.labeled.positive)
    labeled.positive.frac.covered <- sum(replicates.labeled.positive$frac.covered)

    # To estimate precision and recall, we do not know the true
    # positive samples. (Here, true positive for a sample might
    # mean that there's RNA in the region that SHERLOCK is targeting.)
    # But we do know how much of each sample's genome we could
    # sequence, so we can use this. Let's imagine that what we're actually
    # classifying are regions in the Zika genome that SHERLOCK is
    # targeting. An estimate of our true positive count is
    # total.frac.covered, assuming that the distribution of coverage in
    # the Zika genome is uniform. (To see this, note that the estimate
    # of the true positive count is sum_{genome g}(Pr(there is RNA
    # in genome g in the region SHERLOCK is targeting)) = sum_{g}(frac
    # covered in genome g) = total.frac.covered.) The number of these
    # regions that we labeled positive is just num.positive, since we
    # label a sample as positive and each sample has just one such
    # region. The number that we correctly labeled as positive is
    # estimated by labeled.positive.frac.covered (like before, this
    # quantity is estimated by sum_{positive sample p}(Pr(there is
    # RNA in the genome of p in the region SHERLOCK is targeting)) =
    # sum_{p}(frac covered in genome of p) = labeled.positive.frac.covered).

    # Precision = (num correctly labeled positive) / (num labeled positive)
    precision = labeled.positive.frac.covered / num.positive
    # Recall    = (num correctly labeled positive) / (num true positives)
    recall = labeled.positive.frac.covered / total.frac.covered

    F1.score <- 2 * (precision * recall) / (precision + recall)
    F0.5.score <- (1 + 0.5^2) * (precision * recall) / (0.5^2 * precision + recall)

    return(c(min.pool.thres, precision, recall, F1.score, F0.5.score))
}

expected.num.genomes <- function(frac.covered.thres) {
    # Calculate an "expected" number of genomes covered given a
    # threshold, frac.covered.thres, on the fraction of the genome
    # covered
    # Any replicates with frac.covered >= frac.covered.thres are labeled
    # as "positive"
    # Returns (frac.covered.thres, expected num of genomes)

    # Pull out the replicates that are labeled as positive
    replicates.labeled.positive <- replicates[replicates$frac.covered >= frac.covered.thres, ]
    num.positive <- nrow(replicates.labeled.positive)

    expected.num <- num.positive * frac.covered.thres

    return(c(frac.covered.thres, expected.num))
}

########################################
# Find optima of the objective functions
########################################

call.obj.fn <- function(vals.to.test, obj.fn) {
    # Return a dataframe calling obj.fn on vals.to.test
    vals.to.test <- unique(sort(vals.to.test))
    out.vec <- lapply(vals.to.test, obj.fn)
    out.df <- as.data.frame(do.call(rbind, out.vec))
    return(out.df)
}

find.maxima <- function(var.to.solve, obj.var) {
    # Return the values in var.to.solve that correspond to
    # maxima in obj.var; obj.var[i] should correspond to
    # the value of the objective function givne var.to.solve[i]
    return(var.to.solve[which(obj.var == max(obj.var))])
}

# Calculate F1 scores over possible amplicon concentration thresholds
accuracies <- call.obj.fn(replicates$min.pool, accuracy.measure)
colnames(accuracies) <- c("min.pool.thres", "precision", "recall", "F1.score", "F0.5.score")
accuracies$min.pool.thres.plus.one <- accuracies$min.pool.thres + 1
print(paste("Threshold(s) on 1+min(pool concentration) with max F_1 score:",
    find.maxima(accuracies$min.pool.thres.plus.one, accuracies$F1.score)))
print(paste("Threshold(s) on 1+min(pool concentration) with max F_0.5 score:",
    find.maxima(accuracies$min.pool.thres.plus.one, accuracies$F0.5.score)))

# Calculate "expected" num genomes over possible fraction covered thresholds
num.gnm <- call.obj.fn(replicates$frac.covered, expected.num.genomes)
colnames(num.gnm) <- c("frac.covered.thres", "expected.num.genomes")
print(paste("Threshold(s) on frac covered with max expected num genomes:",
    find.maxima(num.gnm$frac.covered.thres, num.gnm$expected.num.genomes)))

############
# Make plots
############

# Plot the objective values versus possible thresholds
p.accuracies <- ggplot(accuracies, aes(x=min.pool.thres.plus.one)) +
    geom_line(aes(y=precision, color="precision"), size=1) +
    geom_line(aes(y=recall, color="recall"), size=1) +
    geom_line(aes(y=F1.score, color="F_1 score"), size=2) +
    geom_line(aes(y=F0.5.score, color="F_0.5 score"), size=2) +
    scale_size(guide="none") +
    scale_x_log10() +
    xlab("Threshold on 1 + min(pool concentration)") +
    ylab("Value") +
    ggtitle("Threshold on pool concentration to maximize F_beta score")

p.num.gnm <- ggplot(num.gnm, aes(x=frac.covered.thres)) +
    geom_line(aes(y=expected.num.genomes), size=1) +
    scale_size(guide="none") +
    xlab("Threshold on fraction of genome covered") +
    ylab("Expected number of genomes") +
    ggtitle("Threshold on genome covered to maximize 'expected' number of genomes")

g <- arrangeGrob(p.accuracies, p.num.gnm)
ggsave(file=out.pdf, g, width=8, height=8, useDingbats=FALSE)

write.table(accuracies, file=out.accuracies.tsv, row.names=FALSE)
write.table(num.gnm, file=out.expected.nnm.gnm.tsv, row.names=FALSE)