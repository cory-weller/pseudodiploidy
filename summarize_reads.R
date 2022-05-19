#!/usr/bin/env Rscript

library(data.table)
library(foreach)

# o <- foreach(i = 1:8, .combine='rbind') %do% {
#     foreach(gene=c('spo13','his2','can1'), .combine='rbind') %do% {
#         filename <- paste0('20220314-S', i,, '-', gene, '.counts')
#         dat <- fread(filename)
#         dat[, sum(V1), by=V2]
#     }
# }