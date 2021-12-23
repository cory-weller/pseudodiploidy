#!/usr/bin/env Rscript

library(data.table)
library(DESeq2)
library(foreach)

args <- commandArgs(trailingOnly = TRUE)

# args <- c("data/input/combined_featurecounts.csv",
#            "data/input/samples.csv",
#            "data/processed/featurecounts.mat.RDS")

featureCountsFilename <- args[1]
coldataFilename <- args[2]
outputFilename <- args[3]

coldata <- fread(coldataFilename)
coldataIDs <- coldata[, id]
coldata[, id := NULL]
coldata <- as.data.frame(coldata)
rownames(coldata) <- coldataIDs

featureCounts <- fread(featureCountsFilename)
featureCounts[, "msp3_6_1_alpha" := featureCounts[,"msp3_6_1"]]
featureCounts[, "msp3_6_2_alpha" := featureCounts[,"msp3_6_2"]]
featureCounts[, "msp3_6_3_alpha" := featureCounts[,"msp3_6_3"]]
featureCounts[, "msp3_6_4_alpha" := featureCounts[,"msp3_6_4"]]
setnames(featureCounts, "msp3_6_1", "msp3_6_1_a")
setnames(featureCounts, "msp3_6_2", "msp3_6_2_a")
setnames(featureCounts, "msp3_6_3", "msp3_6_3_a")
setnames(featureCounts, "msp3_6_4", "msp3_6_4_a")
setcolorder(featureCounts, c("Geneid", "Chr", "Start", "End", "Strand", "Length", rownames(coldata)))
counts <- as.matrix(featureCounts[,!1:6])
rownames(counts) <- featureCounts[,Geneid]

saveRDS(counts, file = outputFilename)

