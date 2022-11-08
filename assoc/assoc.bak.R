#!/usr/bin/env Rscript
library(readr)
library(dplyr)
library(tidyr)

assoc_quant <- function(expr_path, pheno_path) {
    df <- read_delim(expr_path, '\t', col_names=c('PN', 'kmer', 'expr'))
    df <- df %>% inner_join(read_delim(file.path(pheno_path, 'qtl.val'), '\t', col_names=c('PN', 'gt')), by='PN')
    df <- df %>% complete(kmer, nesting(PN, gt), fill=list(expr=0))

    kmers <- pull(arrange(distinct(df, kmer), kmer))
    output <- data.frame(kmer=kmers, p=rep(0, length(kmers)))

    for (km in kmers) {
        ll <- glm(gt ~ expr, family=gaussian(), data=df %>% filter(kmer==km))
        cc <- coef(summary(ll))
        if (length(cc[,'Pr(>|t|)']) > 1) {
            output$p[output$kmer == km] <- cc['expr','Pr(>|t|)']
        }
    }
    return(output[output$p != 0,])
}

assoc_bin <- function(expr_path, pheno_path) {
    df <- bind_rows(
        read_tsv(file.path(pheno_path, 'aff.auto.pnlist'), col_names = c('PN')) %>% mutate(aff=1),
        read_tsv(file.path(pheno_path, 'ctl.auto.pnlist'), col_names = c('PN')) %>% mutate(aff=0),
    ) %>% inner_join(read_tsv(file.path(pheno_path, 'fitted.auto.tsv'), col_names=c('PN','fitted')), by='PN')

    df <- df %>% inner_join(read_delim(expr_path, '\t', col_names = c('PN', 'kmer', 'expr')))
    df <- df %>% complete(kmer, nesting(PN, aff, fitted), fill=list(expr=0))

    kmers <- pull(arrange(distinct(df, kmer), kmer))
    output <- data.frame(kmer=kmers, p=rep(0, length(kmers)))

    for (km in kmers) {
        ll <- glm(aff ~ expr + offset(fitted), family=binomial(), data=df %>% filter(kmer==km))
        cc <- coef(summary(ll))
        if (length(cc[,'Pr(>|z|)']) > 1) {
            output$p[output$kmer == km] <- cc['expr','Pr(>|z|)']
        }
    }
    return(output[output$p != 0,])
}

main <- function(args) {
    if (file.exists(file.path(args[2], 'qtl.val'))) {
        output <- assoc_quant(args[1], args[2])
    } else if (file.exists(file.path(args[2], 'aff.auto.pnlist'))) {
        output <- assoc_bin(args[1], args[2])
    }
    write.table(output, args[3], col.names=F, row.names=F, quote=F)
}

if (sys.nframe() == 0) {
    args <- commandArgs(trailingOnly=T)
    main(args)
}
