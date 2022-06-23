#!/usr/bin/env Rscript

library(data.table)
library(foreach)
library(ggplot2)
library(ggthemes)

o <- foreach(gene=c('his2','spo13','can1'), .combine='rbind', .errorhandling='remove') %do% {
    foreach(sample=c('S1', 'S2', 'S3', 'S4', 'S5', 'S6', 'S7', 'S8'), .combine='rbind', .errorhandling='remove') %do% {
    file <- paste0('data/processed/20220314-', sample, '-', gene, '.counts')
    print(file)
    dt <- fread(file)

    setnames(dt, c('count','type','seq'))
    dt[, 'length' := nchar(seq)]
    out <- dt[, list(sample, gene, 'N'=sum(count)), by=list(type,length)][]
    return(out)
    }
}

o.counts <- o[, c("sample", "gene", "type", "N")]



nomatch <- foreach(gene=c('his2','spo13','can1'), .combine='rbind', .errorhandling='remove') %do% {
    foreach(sample=c('S1', 'S2', 'S3', 'S4', 'S5', 'S6', 'S7', 'S8'), .combine='rbind', .errorhandling='remove') %do% {
    nomatch_file <- paste0('data/processed/20220314-', sample, '-nomatch.counts')
    nomatch_counts <- sum(fread(nomatch_file, select=1)$V1)
    data.table(sample, gene='nomatch', N=nomatch_counts)
    }
}

total_reads <- rbindlist(list(o.counts[,c('sample','gene','N')], nomatch))[, .('reads'=sum(N)), by=sample]

o[sample=='S1', "treatment" := "wt"]
o[sample=='S1', "mating_type" := "a"]
o[sample=='S2', "treatment" := "pseudo"]
o[sample=='S2', "mating_type" := "a"]
o[sample=='S3', "treatment" := "wt"]
o[sample=='S3', "mating_type" := "alpha"]
o[sample=='S4', "treatment" := "pseudo"]
o[sample=='S4', "mating_type" := "alpha"]
o[sample=='S5', "treatment" := "nej1-delete"]
o[sample=='S5', "mating_type" := "a"]
o[sample=='S6', "treatment" := "nej1-delete"]
o[sample=='S6', "mating_type" := "alpha"]
o[sample=='S7', "treatment" := "NAM"]
o[sample=='S7', "mating_type" := "a"]
o[sample=='S8', "treatment" := "NAM"]
o[sample=='S8', "mating_type" := "alpha"]

o[type=='wt', type := 'unedited']

setkey(total_reads, sample)
setkey(o, sample)

o <- merge(o, total_reads)
o[, fraction_reads := N/reads]


 spo13 <- ggplot(o[gene=='spo13' & length < 500], aes(x=length,y=fraction_reads, color=mating_type)) +
  geom_point(shape=21, alpha=0.7) +
  facet_grid(treatment~type) +
    labs(title='spo13', x='read length', y='Fraction of reads (logit-transformed)', color='Mating Type') +
    scale_y_continuous(trans='logit', breaks=c(0,0.2,0.4,0.6)) +
    geom_vline(xintercept=268, alpha=0.5, linetype='dotted') +
    theme_few(12) + guides(color=guide_legend('Mating Type'))



 can1 <- ggplot(o[gene=='can1' & length < 500], aes(x=length,y=fraction_reads, color=mating_type)) +
  geom_point(shape=21, alpha=0.7) +
  facet_grid(treatment~type) +
    labs(title='can1', x='read length', y='Fraction of reads (logit-transformed)') +
    scale_y_continuous(trans='logit', breaks=c(0,0.2,0.4,0.6)) +
    geom_vline(xintercept=315, alpha=0.5, linetype='dotted') +
    theme_few(12)




#  his2 <- ggplot(o[gene=='his2' & length < 500], aes(x=length,y=fraction_reads, color=mating_type)) +
#   geom_point(shape=21, alpha=0.7) +
#   facet_grid(treatment~type) +
#     labs(title='his2', x='read length', y='Fraction of reads (logit-transformed)') +
#     scale_y_continuous(trans='logit', breaks=c(0,0.2,0.4,0.6)) +
#     geom_vline(xintercept=280, alpha=0.5, linetype='dotted') +
#     theme_few(14)

  ggsave(spo13, file='spo13.png', width=24, height=24, units='cm')
  ggsave(can1, file='can1.png', width=24, height=24, units='cm')
  # ggsave(his2, file='his2.png', width=24, height=12, units='cm')
