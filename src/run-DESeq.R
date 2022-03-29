#!/usr/bin/env Rscript

library(data.table)
library(DESeq2)
library(foreach)

args <- commandArgs(trailingOnly = TRUE)


featureCountsFilename <- args[1]
coldataFilename <- args[2]
outputFilename <- args[3]

coldata <- fread(coldataFilename)
coldataIDs <- coldata[, id]
coldata[, id := NULL]
coldata <- as.data.frame(coldata)
rownames(coldata) <- coldataIDs


dds <- DESeqDataSetFromMatrix(countData=counts, colData=coldata, design = ~background*treatment)

# Filter out rows (genes) with counts fewer than 10 across all samples
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

print("Running DESeq")
dds <- DESeq(dds)

saveRDS(dds, file = outputFilename)


dds$group <- factor(paste0(dds$background, dds$treatment))
design(dds) <- ~ group
dds <- DESeq(dds)

#
# > resultsNames(dds)
# [1] "Intercept"                       "group_alphadiploid_vs_adiploid"
# [3] "group_alphanam_vs_adiploid"      "group_alphanej1_del_vs_adiploid"
# [5] "group_alphano_vs_adiploid"       "group_alphapseudo_vs_adiploid"
# [7] "group_anam_vs_adiploid"          "group_anej1_del_vs_adiploid"
# [9] "group_ano_vs_adiploid"           "group_apseudo_vs_adiploid"
# 



## Looking at TPM
dat <- fread(featureCountsFilename)
TPM.mat <- foreach(sample=colnames(dat)[7:length(colnames(dat))], .combine="cbind") %do% {
    counts <- dat[[sample]]
    rates <- counts / dat[,Length]
    TPM <- data.table((rates / sum(rates)) * 1e6)
    setnames(TPM, sample)
    return(TPM)
}
TPM.mat [, Geneid := dat[,Geneid]]
TPM.long <- melt(TPM.mat, measure.vars=colnames(dat)[7:length(colnames(dat))], variable.name="id")
coldata <- fread(coldataFilename)
setkey(coldata, id)
setkey(TPM.long, id)

TPM.long <- merge(TPM.long, coldata, all=TRUE)
TPM.long[, type := NULL]
TPM.long[is.na(treatment), background := "diploid"]
TPM.long[is.na(treatment), treatment := "diploid"]
# Plot expression of NEJ1 across samples
ggplot(data=TPM.long[Geneid == "YLR265C"], mapping=aes(x=treatment, y=value, color=background)) +
geom_point(shape=21, alpha=0.8, size=3) +
labs(x="Treatment", y="Transcripts per Million (TPM)", color="Background") +
theme_few(14)
fwrite(TPM.long, file="TPM.long.txt", quote=F, row.names=F, col.names=T, sep="\t")

setcolorder(TPM.mat, c("Geneid", colnames(TPM.mat))[1:(length(colnames(TPM.mat))-1)])
fwrite(TPM.mat, file="TPM.txt", quote=F, row.names=F, col.names=T, sep="\t")

#resultsNames(dds)
#res <- results(dds, name="")

runContrast <- function(group_factor, group1, group2) {
    res.1 <- results(dds, contrast = c(group_factor, group1, group2))
    res <- lfcShrink(dds, contrast=c(group_factor, group1, group2), res=res.1)
    res <- setDT(as.data.frame(res), keep.rownames=T)

    # plot

    library(ggplot2)
    library(ggrepel)
    library(ggthemes)

    res[, plotAlpha := ifelse(abs(log2FoldChange) > 1, 1, 0.5)]
    res[plotAlpha == 1, labelText := rn]

    g <-  ggplot(res, aes(
                    x=log2FoldChange, 
                    y=-1*log10(padj), 
                    alpha=plotAlpha,
                    label=labelText
                    )) +
    geom_point() +
    geom_text_repel() +
    labs(y="-log10(p)", title=paste(group_factor, group1, group2, sep="_")) +
    guides(alpha="none") +
    theme_few()

    genesUp <- res[log2FoldChange > 1 & padj < 0.0001]
    genesDown <- res[log2FoldChange < -1 & padj < 0.0001]
    return(list(g, genesUp, genesDown))
}

# get features file from SGD
system("")
SGDFeatures <- fread("data/external/SGD_features.tab")
SGDFeatures <- SGDFeatures[,c("V4","V1")]
setnames(SGDFeatures, c("SystematicID", "SGD_ID"))
setkey(SGDFeatures, SystematicID)

nam.v.pseudo <- runContrast("group", "anam", "apseudo")

write(paste("SGD", SGDFeatures[nam.v.pseudo[[2]], SGD_ID], sep=":"), file="enriched-in-nam-vs-pseudo.txt")
write(paste("SGD", SGDFeatures[nam.v.pseudo[[3]], SGD_ID], sep=":"), file="enriched-in-pseudo-vs-nam.txt")



pseudo.v.no <- runContrast("group", "apseudo", "ano")

write(paste("SGD", SGDFeatures[pseudo.v.no[[2]], SGD_ID], sep=":"), file="enriched-in-pseudo-vs-no.txt")
write(paste("SGD", SGDFeatures[pseudo.v.no[[3]], SGD_ID], sep=":"), file="enriched-in-no-vs-pseudo.txt")



write(SGDFeatures[nam.v.pseudo[[2]], SGD_ID], file="enriched-in-nam-vs-pseudo.txt")
write(SGDFeatures[nam.v.pseudo[[3]], SGD_ID], file="enriched-in-pseudo-vs-nam.txt")

nam.v.pseudo[[3]]


nam.v.no <- plotContrast("group", "anam", "ano")
pseudo.v.no <- plotContrast("group", "apseudo", "ano")


