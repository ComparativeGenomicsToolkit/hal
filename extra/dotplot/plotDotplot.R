#!/usr/bin/env Rscript

args <- commandArgs(TRUE)
if (length(args) < 2) {
   print("Usage: plotDotplot.R inputDotplotTsv outputPdf")
   quit(status = 1)
}
inputFile <- args[[1]]
outputPdfFile <- args[[2]]
dots <- read.table(inputFile, header=T)
pdf(outputPdfFile)
print(plot(dots, type="l"))
dev.off()