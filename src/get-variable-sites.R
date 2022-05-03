#!/usr/bin/env Rscript

library(data.table)

args <- commandArgs(trailingOnly=TRUE)
# args <- c('data/external/chromosome12.vcf.gz','data/processed/pooled-distinct-strains.txt','12','171320','172320')
vcfFile <- args[1]
strainsFile <- args[2]
chromosome <- paste0('chromosome', args[3])
start <- as.numeric(args[4])
end <- as.numeric(args[5])


# Load VCF
vcf <- fread(vcfFile, header=FALSE)
headerText <- colnames(fread('data/input/vcf-header.txt'))
setnames(vcf, headerText)
setnames(vcf, '#CHROM', 'CHROM')

strains <- fread(strainsFile, header=FALSE)$V1
strains <- strains[strains != 'S288C']
nStrains <- length(strains)

desiredCols <- c('CHROM', 'POS', 'REF', 'ALT', strains)
vcf <- vcf[, desiredCols, with=F][POS >= start & POS <= end & CHROM == chromosome]
vcf[, nRef := apply(.SD, 1, function(x) sum(x == '0/0')), .SDcols=strains]
vcf[, nAlt1 := apply(.SD, 1, function(x) sum(x == '1/1')), .SDcols=strains]
vcf[, nAlt2 := apply(.SD, 1, function(x) sum(x == '2/2')), .SDcols=strains]
vcf[, nAlt3 := apply(.SD, 1, function(x) sum(x == '3/3')), .SDcols=strains]
vcf[, nAlt4 := apply(.SD, 1, function(x) sum(x == '4/4')), .SDcols=strains]
vcf[, nAlt5 := apply(.SD, 1, function(x) sum(x == '5/5')), .SDcols=strains]
vcf[, nAlt6 := apply(.SD, 1, function(x) sum(x == '6/6')), .SDcols=strains]
vcf[, nAlt7 := apply(.SD, 1, function(x) sum(x == '7/7')), .SDcols=strains]
vcf[, nAlt8 := apply(.SD, 1, function(x) sum(x == '8/8')), .SDcols=strains]
vcf[, nMissing := apply(.SD, 1, function(x) sum(x == './.')), .SDcols=strains]


# exclude all-ref strains
vcf <- vcf[nRef + nMissing != nStrains]
# print variable sites
cat(vcf$POS, sep="\n")