#!/usr/bin/env Rscript

library(data.table)
library(ggplot2)
library(ggthemes)
library(viridis)
library(R.utils)

dat <- fread("data/processed/TPM.txt.gz")

geneOrder <- dat[, sum(value), by=Geneid][order(V1)][,Geneid]

dat[, Geneid := factor(Geneid, levels=geneOrder)]

g <- ggplot(data=dat,
        aes(
            y=treatment,
            x=Geneid,
            fill=log2(value)
        )
    ) +
    geom_tile() +
    facet_grid(background~.) +
    scale_fill_viridis() +
    labs(
            y="Treatment",
            x="GeneID",
            fill="Log2(TPM)"
    ) +
    theme_few() +
    theme(
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()
        )

ggsave(g, file="reports/heatmap.png")