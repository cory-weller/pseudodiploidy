#!/usr/bin/env Rscript
library(data.table)
library(foreach)

args <- commandArgs(trailingOnly = TRUE)


featureCountsFilename <- args[1]
coldataFilename <- args[2]
outputFilename <- args[3]

dat <- fread(featureCountsFilename)
TPM.mat <- foreach(sample=colnames(dat)[7:length(colnames(dat))], .combine="cbind") %do% {
    counts <- dat[[sample]]
    rates <- counts / dat[,Length]
    TPM <- data.table((rates / sum(rates)) * 1e6)
    setnames(TPM, sample)
    return(TPM)
}

TPM.mat[, Geneid := dat[,Geneid]]
TPM.long <- melt(TPM.mat, measure.vars=colnames(dat)[7:length(colnames(dat))], variable.name="id")
coldata <- fread(coldataFilename)
setkey(coldata, id)
setkey(TPM.long, id)

TPM.long <- merge(TPM.long, coldata, all=TRUE)
TPM.long[, type := NULL]
TPM.long[is.na(treatment), background := "diploid"]
TPM.long[is.na(treatment), treatment := "diploid"]
fwrite(TPM.long, file = outputFilename, quote=F, row.names=F, col.names=T, sep="\t")