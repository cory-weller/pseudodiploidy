#!/usr/bin/env Rscript

library(data.table)
library(DESeq2)
library(ggplot2)
library(ggthemes)
library(R.utils)

source('src/vars.R')


runContrast <- function(dds, group1, group2) {
    if(levels(dds$condition)[1] != group1) {
        dds$condition <- relevel(dds$condition, ref=group1)
        dds <- estimateSizeFactors(dds)
        dds <- estimateDispersions(dds)
        dds <- nbinomWaldTest(dds)
    }
    res <- lfcShrink(  dds = dds,
                coef = paste('condition_', group2, '_vs_', group1, sep=''),
                type = 'apeglm'
    )
    return(as.data.table(res))
}

# Load in DESeq2 data set (DDS) object
dds <- readRDS(DESeqDdsFilename)

nam.v.pseudo <- runContrast(dds, "a_nam", "a_pseudo")


# get features file from SGD
SgdFeatures <- fread(SgdFeaturesFilename, header = FALSE)
SgdFeatures <- SgdFeatures[,c("V4","V1")]
setnames(SgdFeatures, c("Systematic_ID", "SGD_ID"))
setkey(SgdFeatures, Systematic_ID)

