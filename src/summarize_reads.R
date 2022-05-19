#!/usr/bin/env Rscript

library(data.table)
library(foreach)

o <- foreach(i = 1:8, .combine='rbind') %do% {
    foreach(gene=c('spo13','his2','can1'), .combine='rbind') %do% {
        filename <- paste0('data/processed/20220314-S', i, '-', gene, '.counts')
        dat <- fread(filename, header=FALSE)
        if(nrow(dat) == 0) {
            return(NULL)
        } else {
            setnames(dat, "V1", "N")
            setnames(dat, "V2", "class")
            dat[, sample := i]
            dat[, gene := gene]
            dat.ag <- dat[, list("N"=sum(N)), by=list(sample, gene, class)]
            return(dat.ag[])
        }
    }
}

library(ggplot2)


o[sample==1, "treatment" := "wt"]
o[sample==1, "mating_type" := "a"]

o[sample==2, "treatment" := "pseudo"]
o[sample==2, "mating_type" := "a"]

o[sample==3, "treatment" := "wt"]
o[sample==3, "mating_type" := "alpha"]

o[sample==4, "treatment" := "pseudo"]
o[sample==4, "mating_type" := "alpha"]

o[sample==5, "treatment" := "nej1-delete"]
o[sample==5, "mating_type" := "a"]

o[sample==6, "treatment" := "nej1-delete"]
o[sample==6, "mating_type" := "alpha"]

o[sample==7, "treatment" := "NAM"]
o[sample==7, "mating_type" := "a"]

o[sample==8, "treatment" := "NAM"]
o[sample==8, "mating_type" := "alpha"]

o[, treatment := factor(treatment, levels=c('wt','nej1-delete','pseudo','NAM'))]
o[, mating_type := factor(mating_type, levels=c('a','alpha'))]
o[, pct := N/sum(N), by=list(gene,mating_type,treatment)]

g <- ggplot(o[gene != 'his2'], aes(x=gene, y=N, fill=class)) + geom_bar(stat='identity',position='dodge') + facet_grid(mating_type~treatment)


ggplot(o[gene != 'his2'], aes(x=gene, y=pct, fill=class)) + geom_bar(stat='identity') + facet_grid(mating_type~treatment)


ggplot(stacked, aes(x=gene, y=percent, fill=