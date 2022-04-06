#!/usr/bin/env R

library(data.table)

guides <- fread('data/processed/chr12-guides.tsv')
setnames(guides, c("guideSeq", "PAM", "start", "stop", "strand"))
variableSites <- fread('data/processed/S288C-chr12-variable-sites.txt')$V1

putativeGuides <- guides[start %between% c(171320, 172320)]