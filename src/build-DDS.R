#!/usr/bin/env Rscript

library(data.table)
library(ggplot2)
library(ggthemes)
library(foreach)
library(R.utils)
library(DESeq2)

source('src/vars.R')

counts <- readRDS(featureCountMatrixFilename)

coldata <- fread(coldataFilename)
coldataIDs <- coldata[, id]
coldata[, id := NULL]
coldata[, condition := paste(background, treatment, sep="_")]
coldata[, c("background", "treatment", "type") := NULL]
coldata <- as.data.frame(coldata)
rownames(coldata) <- coldataIDs


dds <- DESeqDataSetFromMatrix(countData=counts, colData=coldata, design = ~condition)
dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds)
dds <- nbinomWaldTest(dds)

saveRDS(dds, DESeqDdsFilename)