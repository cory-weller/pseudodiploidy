#!/usr/bin/env Rscript
options(width=900)
args <- commandArgs(trailingOnly=TRUE)
#args <- c('data/processed/chr12-guides.tsv', 'data/processed/pooled-distinct-strains-variable-sites.txt', '171320', '172320')
library(data.table)

allGuides <- fread(args[1], header=F)
setnames(allGuides, c('seq', 'PAM', 'startPosition', 'endPosition', 'strand'))
allGuides[, idx := 1:.N]
variableSites <- fread(args[2])$V1
start <- as.numeric(args[3])
end <- as.numeric(args[4])

usableGuides <- allGuides[startPosition >= start & endPosition <= end]

usableGuides[, varSitesBetween := sum(startPosition <= variableSites & endPosition >= variableSites), by=idx]
print(usableGuides)