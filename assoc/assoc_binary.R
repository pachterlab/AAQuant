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

bootstrap <- function(fn, N, df) {

    vec <- c()

    for (i in 1:N) {

        pval <- tryCatch(
                    { fn(expr ~ aff, data=df[sample(1:nrow(df), size=nrow(df), replace=T),])['p.value'] },
                    error=function(cond) {
                        NA
                    },
                    warning=function(cond) {
                        NA
                    }
                )
        vec <- c(vec, as.numeric(pval))

    }
    return(mean(vec, na.rm=T))
}

assoc_bin <- function(expr_path, aff_path, ctl_path) {
    print("assoc_bin")
    df <- bind_rows(
        read_tsv(aff_path, col_names = c('PN')) %>% mutate(aff=1),
        read_tsv(ctl_path, col_names = c('PN')) %>% mutate(aff=0),
    )# %>% inner_join(read_tsv(file.path(pheno_path, 'fitted.auto.tsv'), col_names=c('PN','fitted')), by='PN')

    df$aff <- as.factor(df$aff)

    df <- df %>% inner_join(read_delim(expr_path, '\t', col_names = c('PN', 'kmer', 'expr')))
    kmers <- pull(arrange(distinct(df, kmer), kmer))
    output <- data.frame(kmer=kmers, p_t_test=rep(0, length(kmers)), p_wilcoxon=rep(0, length(kmers)))
    if (nrow(df) == 0) {
        return(output)
    }

    df <- df %>% complete(kmer, nesting(PN, aff), fill=list(expr=0))


    for (km in kmers) {
        # ll <- t.test(expr ~ aff, data=df %>% filter(kmer==km))
        # output$p_t_test[output$kmer == km] <- as.numeric(ll['p.value'])
        ll <- wilcox.test(expr ~ aff, data=df %>% filter(kmer==km))
        output$p_wilcoxon[output$kmer == km] <- as.numeric(ll['p.value'])

        # pval <- tryCatch(
        #             { wilcox.test(expr ~ aff, data=df %>% filter(kmer==km))['p.value'] },
        #             error=function(cond) {
        #                 1
        #             },
        #             warning=function(cond) {
        #                 1
        #             }
        #         )
        # output$p_wilcoxon <- as.numeric(pval)

        # output$p_t_test[output$kmer == km] <- bootstrap(t.test, 100, df %>% filter(kmer==km))
        # output$p_wilcoxon[output$kmer == km] <- bootstrap(wilcox.test, 100, df %>% filter(kmer==km))
    }
    return(output[output$p_wilcoxon != 0,])
}

main <- function(args) {
    output <- assoc_bin(args[1], args[2], args[3])
    write.table(output, args[4], col.names=F, row.names=F, quote=F)
}

if (sys.nframe() == 0) {
    args <- commandArgs(trailingOnly=T)
    print(args)
    main(args)
}
