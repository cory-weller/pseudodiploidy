#!/usr/bin/env Rscript

library(data.table)
library(DESeq2)
library(ggplot2)
library(ggthemes)
library(R.utils)

source('src/vars.R')

args <- commandArgs(trailingOnly = TRUE)
group1 <- args[1]
group2 <- args[2]

contrasts <- fread(contrastsFilename, header = FALSE)
nContrasts <- nrow(contrasts)
bonferroni <- 0.05 / nContrasts

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

DEgenes <- data.table()

# Load in DESeq2 data set (DDS) object
dds <- readRDS(DESeqDdsFilename)

for(i in 1:nContrasts) {
    group1 <- contrasts[i,V1]
    group2 <- contrasts[i,V2]

    expressionContrast <- runContrast(dds, group1, group2)
    # correct contrast direction
    expressionContrast[, log2FoldChange := -1*log2FoldChange]

    # add isSignif column, where TRUE if log2FoldChange and padj pass threshold
    expressionContrast[, isSignif := FALSE]
    expressionContrast[, plotColor := 'black']
    #expressionContrast[, geneLabel := '']
    expressionContrast[, plotAlpha := 0.3]
    expressionContrast[abs(log2FoldChange) > foldChangeSignificance & padj < bonferroni, isSignif := TRUE]

    # add gene label column only for significant genes
    #expressionContrast[isSignif==TRUE, geneLabel := geneID]
    expressionContrast[isSignif==TRUE, plotAlpha := 1.0]
    expressionContrast[isSignif==TRUE, plotColor := 'red']
    signifDEgenes <- expressionContrast[isSignif==TRUE]
    signifDEgenes[, 'group1' := group1]
    signifDEgenes[, 'group2' := group2]
    
    DEgenes <- rbindlist(list(DEgenes, signifDEgenes))
    backgroundStrain <- unlist(strsplit(group1, '_'))[1]

    g <- ggplot(data=expressionContrast, 
            aes(
                x=log2FoldChange, 
                y=-1*log10(padj), 
                alpha=plotAlpha
            )
        ) + 
        geom_point(shape=21, color=expressionContrast$plotColor) +
        labs(x=paste('log2 fold change in ', group1, ' relative to haploid ', backgroundStrain, sep='')) +
        theme_few(12) +
        guides(alpha='none') +
        guides(color='none') +
        geom_hline(yintercept = -1*log10(bonferroni), linetype='dashed', alpha=0.5, color='gray') +
        geom_vline(xintercept = foldChangeSignificance, linetype='dashed', alpha=0.5, color='gray') +
        geom_vline(xintercept = -foldChangeSignificance, linetype='dashed', alpha=0.5, color='gray')

    ggsave(g, file=paste0(group1, '-vs-', group2, '.png'), width=20, height=30, units='cm')
}
DEgenes[, minusLog10p := -1*log10(padj)]
DEgenes[, c('plotAlpha', 'plotColor', 'isSignif', 'padj') := NULL]

fwrite(DEgenes, file=DEgenesFilename, quote=F, row.names=F, col.names=T, sep='\t')

# get list of upregulated genes
# expressionContrast[log2FoldChange > 2 & padj < 0.05, geneID]

# get list of downregulated genes
# expressionContrast[log2FoldChange < -2 & padj < 0.05, geneID]



# TODO: Code below to merge in SGD features
# get features file from SGD
# SgdFeatures <- fread(SgdFeaturesFilename, header = FALSE)
# SgdFeatures <- SgdFeatures[,c("V4","V1")]
# setnames(SgdFeatures, c("Systematic_ID", "SGD_ID"))
# setkey(SgdFeatures, Systematic_ID)