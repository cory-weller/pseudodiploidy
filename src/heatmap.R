#!/usr/bin/env Rscript

library(data.table)
library(ggplot2)
library(ggthemes)
library(viridis)
library(R.utils)

dat <- fread("data/processed/TPM.txt.gz")

geneOrder <- dat[, sum(value), by=Geneid][order(V1)][,Geneid]

dat[, Geneid := factor(Geneid, levels=geneOrder)]

# average TPM across treatments
dat.ag <- dat[, list("meanTPM" = mean(value, na.rm=T)), by=list(Geneid, background, treatment)]
dat.dip <- dat.ag[background=="diploid"]
dat.dip[, c("background", "treatment") := NULL]
setnames(dat.dip, 'meanTPM', 'diploidTPM')
setkey(dat.dip, Geneid)
setkey(dat.ag, Geneid)

dat.merge <- merge(dat.ag, dat.dip)

dat.merge[, FC_rel_to_diploid := meanTPM / diploidTPM]


g <- ggplot(data=dat.merge[background != 'diploid'],
        aes(
            y=background,
            x=Geneid,
            fill=log2(FC_rel_to_diploid)
        )
    ) +
    geom_tile() +
    facet_grid(.~treatment) +
    labs(
            y="Treatment",
            x="GeneID",
            fill="Log2(TPM) relative to diploid"
    ) +
    theme_few() +
    theme(
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()
        ) +
    coord_flip() +
    scale_fill_gradientn(colours = c("cyan", "black", "red"),
                       values = scales::rescale(c(-10, -.5, 0, .5, 14)))

ggsave(g, file="reports/heatmap.png")