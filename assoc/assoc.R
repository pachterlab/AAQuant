#!/usr/bin/env Rscript

#.libPaths("/home/kristjanh/R/x86_64-redhat-linux-gnu-library/3.6")
#print(.libPaths())
#print(R.version.string)
library(readr)
library(dplyr)
library(tidyr)
library(harmonicmeanp)

rank_transform <- function(x) {
    return(qnorm((rank(x,ties.method='random')-0.5)/length(x)))
}

assoc_quant <- function(expr_path, pheno_path) {
    print("assoc_quant")
    df_pheno <- read_delim(pheno_path, col_names=c('PN', 'gt'))
    df_pheno$gt <- rank_transform(df_pheno$gt)
    df <- read_delim(expr_path, '\t', col_names=c('PN', 'kmer', 'expr'))
    df <- df %>% inner_join(df_pheno, by='PN')
    df <- df %>% complete(kmer, nesting(PN, gt), fill=list(expr=0))

    kmers <- pull(arrange(distinct(df, kmer), kmer))
    output <- data.frame(kmer=kmers, p=rep(0, length(kmers)))

    for (km in kmers) {
        ll <- glm(gt ~ expr-1, family=gaussian(), data=df %>% filter(kmer==km))
        cc <- coef(summary(ll))
        if (length(cc[,'Pr(>|t|)']) >= 1) {
            output$p[output$kmer == km] <- cc['expr','Pr(>|t|)']
        }
    }
    return(output[output$p != 0,])
}

assoc_quant_add <- function(expr_path, pheno_path) {
    print("assoc_quant_add")
    df_pheno <- read_delim(pheno_path, col_names=c('PN', 'gt'))
    df <- read_delim(expr_path, '\t', col_names=c('PN', 'kmer', 'expr'))
    df <- df %>% inner_join(df_pheno, by='PN')
    df <- df %>% complete(kmer, nesting(PN, gt), fill=list(expr=0))

    kmers <- pull(arrange(distinct(df, kmer), kmer))
    output <- data.frame(kmer=kmers, p=rep(0, length(kmers)))

    for (km in kmers) {
        ll <- glm(gt ~ expr, family=gaussian(), data=df %>% filter(kmer==km))
        cc <- coef(summary(ll))
        if (length(cc[,'Pr(>|t|)']) >= 1) {
            output$p[output$kmer == km] <- cc['expr','Pr(>|t|)']
        }
    }
    return(output[output$p != 0,])
}

assoc_bin <- function(expr_path, pheno_path) {
    print("assoc_bin")
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

boostrap_quant_assoc <- function(df, N) {

    vec <- c()

    for (i in 1:N) {

        ll <- glm(gt ~ expr, family=gaussian(), data=df[sample(1:nrow(df), size=nrow(df), replace=T),])
        cc <- coef(summary(ll))

        if (length(cc[,'Pr(>|t|)']) > 1) {
            if (cc['expr', 'Pr(>|t|)'] > 0) {
                vec <- c(vec, cc['expr','Pr(>|t|)'])
            }
        }

    }
    pval <- tryCatch(
                { p.hmp(vec, w=rep(1/length(vec), length(vec)), L=length(vec)) },
                error=function(cond) {
                    return(min(vec))
                },
                warning=function(cond) {
                    return(min(vec))
                }
            )

    return(pval)
    # return(mean(vec, na.rm=T))

}

boostrap_quant_assoc_driver <- function(expr_path, pheno_path) {
    print("bootstrap_quant_assoc_driver")
    df_pheno <- read_delim(pheno_path, col_names=c('PN', 'gt'))
    df_pheno$gt <- df_pheno$gt * 0.001
    # df_pheno$gt <- rank_transform(df_pheno$gt)
    df <- read_delim(expr_path, '\t', col_names=c('PN', 'kmer', 'expr'))
    df <- df %>% inner_join(df_pheno, by='PN')
    df <- df %>% complete(kmer, nesting(PN, gt), fill=list(expr=0))

    kmers <- pull(arrange(distinct(df, kmer), kmer))
    output <- data.frame(kmer=kmers, p=rep(0, length(kmers)))

    for (km in kmers) {
        output$p[output$kmer == km] <- boostrap_quant_assoc(df %>% filter(kmer==km), 100)
    }
    return(output[output$p != 0,])
}

main <- function(args) {
    if (file.exists(file.path(args[2], 'qtl.val'))) {
        output <- assoc_quant(args[1], file.path(args[2], 'qtl.val'))
    } else if (file.exists(file.path(args[2], 'work/qtl.raw'))) {
        output <- assoc_quant(args[1], file.path(args[2], 'work/qtl.raw'))
    } else if (file.exists(file.path(args[2], 'aff.auto.pnlist'))) {
        output <- assoc_bin(args[1], args[2])
    } else {
        # Run associations against home-made phenotypes. They are assumed to 
        # be quantitative.
        # output <- assoc_quant_add(args[1], args[2])
        # output <- assoc_quant(args[1], args[2])
        output <- boostrap_quant_assoc_driver(args[1], args[2])
    }
    write.table(output, args[3], col.names=F, row.names=F, quote=F)
}

if (sys.nframe() == 0) {
    args <- commandArgs(trailingOnly=T)
    print(args)
    main(args)
}
