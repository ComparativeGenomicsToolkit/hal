#!/usr/bin/env Rscript
# Usage: plotInsertionStats.R <insertionStats.tsv>
# Produces: filename.pdf: density of insertion size in each genome
#           filename-ecdf.pdf: cumulative % of insertions at or below X insertion size in each genome
#           filename-seqFrac-density.pdf: density plot of insertion size, weighted by insertion size
#           filename-seqFrac-hist.pdf: histogram of insertion size, weighted by insertion size (i.e. height of bar = fraction of total inserted sequence)
require(ggplot2)
require(data.table)

args <- commandArgs(TRUE)
insertionStatsFile <- args[[1]]

#ancestorNames <- list(Anc8="human-chimp", Anc7="human-rhesus", Anc5="euarchontoglires", Anc6="mouse-rat", Anc0="boreoeutherian", Anc1="laurasiatheria", Anc2="pig-cow-horse", Anc3="pig-cow", Anc4="dog-cat")

insertionStats <- read.table("entirelyInsertedSequences.tsv", sep="\t", header=T)
insertionStats$maskedBasesFraction <- insertionStats$maskedBases / insertionStats$insertionSize

# Rename the 10-way ancestors
## insertionStats$genome = sapply(insertionStats$genome, function(x) {
##     if (x %in% names(ancestorNames)) {
##         ancestorNames[[as.character(x)]]
##     } else {
##         as.character(x)
##     }})

plot <- ggplot(insertionStats, aes(x=insertionSize)) + scale_x_log10() + facet_wrap(~ genome) + theme_bw()

pdf(paste(insertionStatsFile, ".pdf", sep=""), width=10, height=10)
print(plot + geom_density() + ylab("Density"))
dev.off()

pdf(paste(insertionStatsFile, "-ecdf.pdf", sep=""), width=10, height=10)
print(plot + stat_ecdf() + ylab("Cumulative fraction of insertions"))
dev.off()

pdf(paste(insertionStatsFile, "-seqFrac-density.pdf", sep=""), width=10, height=10)
# stolen from stackoverflow -- get per-genome total insertion size so we can weight correctly.
insertionStats.dt <- data.table(insertionStats)
insertionStats.dt[, totalInsertionSize.per.genome := sum(insertionSize), genome]
print(ggplot(insertionStats.dt, aes(x=insertionSize, weight=insertionSize/totalInsertionSize.per.genome)) + scale_x_log10() + facet_wrap(~ genome) + theme_bw() + geom_density() + ylab("Density weighted by insertion size"))
dev.off()

pdf(paste(insertionStatsFile, "-seqFrac-hist.pdf", sep=""), width=10, height=10)
print(ggplot(insertionStats.dt, aes(x=insertionSize, weight=insertionSize/totalInsertionSize.per.genome)) + scale_x_log10() + facet_wrap(~ genome) + theme_bw() + geom_histogram() + ylab("Fraction of total inserted bases in genome"))
dev.off()
