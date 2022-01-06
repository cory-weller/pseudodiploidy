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
    return(as.data.table(res, keep.rownames="geneID"))
}

# Load in DESeq2 data set (DDS) object
dds <- readRDS(DESeqDdsFilename)

group1 <- "a_nam"
group2 <- "a_pseudo"

expressionContrast <- runContrast(dds, group1, group2)

# add isSignif column, where TRUE if log2FoldChange and padj pass threshold
expressionContrast[abs(log2FoldChange) > 2 & padj < 0.05, isSignif := TRUE]

# add gene label column only for significant genes
expressionContrast[isSignif==TRUE, geneLabel := geneID]

# make column for scatter plot alpha, where insignif genes are see-through
expressionContrast[, plotAlpha := 0.3]
expressionContrast[isSignif==TRUE, plotAlpha := 1.0]

ggplot(data=expressionContrast, 
        aes(
            x=log2FoldChange, 
            y=-1*log10(padj), 
            label=geneLabel, 
            alpha=plotAlpha
        ) + 
        geom_point(shape=21, alpha=0.3) +
        geom_text()


# get list of upregulated genes
expressionContrast[log2FoldChange > 2 & padj < 0.05, geneID]

# get list of downregulated genes
expressionContrast[log2FoldChange < -2 & padj < 0.05, geneID]



# TODO: Code below to merge in SGD features
# get features file from SGD
SgdFeatures <- fread(SgdFeaturesFilename, header = FALSE)
SgdFeatures <- SgdFeatures[,c("V4","V1")]
setnames(SgdFeatures, c("Systematic_ID", "SGD_ID"))
setkey(SgdFeatures, Systematic_ID)

