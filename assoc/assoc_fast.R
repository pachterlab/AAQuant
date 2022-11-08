#!/nfs/fs1/bioinfo/apps-x86_64/R/3.6.1/bin/Rscript
print(.libPaths())
print(R.version.string)
library(readr)
library(dplyr)
library(tidyr)
library(RcppArmadillo)

rank_transform <- function(x){
    return(qnorm((rank(x,ties.method='random')-0.5)/length(x)))
}

assoc_fast <- function(expr_path, pheno_path) {
    print("assoc_fast")
    df_pheno <- read_delim(pheno_path, '\t', col_names=c('PN', 'gt'))
    df <- read_delim(expr_path, '\t', col_names=c('PN', 'kmer', 'expr'))
    df <- df %>% inner_join(df_pheno, by='PN')
    df <- df %>% complete(kmer, nesting(PN, gt), fill=list(expr=0))
    kmers <- pull(arrange(distinct(df, kmer), kmer))
    output <- data.frame(kmer=kmers, coef=rep(0, length(kmers)), p=rep(0, length(kmers)))

    for (km in kmers) {
        dx <- filter(df, kmer==km)
        mdl <- RcppArmadillo::fastLmPure(cbind(1,dx$expr),dx$gt)
        output$coef[output$kmer == km] <- mdl$coefficients[2,1]
        output$p[output$kmer == km] <- 2*pt(abs(mdl$coefficients/mdl$stderr), mdl$df.residual, lower.tail=F)[2,1]
    }
    return(output[output$p != 0,])
}

assoc_fast_rank_transform <- function(expr_path, pheno_path) {
    print("assoc_fast_rank_transform")
    df_pheno <- read_delim(pheno_path, '\t', col_names=c('PN', 'gt'))
    df_pheno$gt <- rank_transform(df_pheno$gt)
    df <- read_delim(expr_path, '\t', col_names=c('PN', 'kmer', 'expr'))
    kmers <- pull(arrange(distinct(df, kmer), kmer))
    df <- df %>% inner_join(df_pheno, by='PN')
    df <- df %>% complete(kmer, nesting(PN, gt), fill=list(expr=0))
    output <- data.frame(kmer=kmers, coef=rep(0, length(kmers)), p=rep(0, length(kmers)))

    for (km in kmers) {
        dx <- filter(df, kmer==km)
        mdl <- RcppArmadillo::fastLmPure(cbind(0, dx$expr),dx$gt)
        output$coef[output$kmer == km] <- mdl$coefficients[2,1]
        output$p[output$kmer == km] <- 2*pt(abs(mdl$coefficients/mdl$stderr), mdl$df.residual, lower.tail=F)[2,1]
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
        ll <- glm(aff ~ expr + offset(fitted) - 1, family=binomial(), data=df %>% filter(kmer==km))
        cc <- coef(summary(ll))
        if (length(cc[,'Pr(>|z|)']) > 1) {
            output$p[output$kmer == km] <- cc['expr','Pr(>|z|)']
        }
    }
    return(output[output$p != 0,])
}

main <- function(args) {
    if (file.exists(file.path(args[2], 'qtl.val'))) {
        output <- assoc_fast_rank_transform(args[1], file.path(args[2], 'qtl.val'))
    } else if (file.exists(file.path(args[2], 'work/qtl.raw'))) {
        output <- assoc_fast_rank_transform(args[1], file.path(args[2], 'work/qtl.raw'))
    } else if (file.exists(file.path(args[2], 'aff.auto.pnlist'))) {
        output <- assoc_bin(args[1], args[2])
    } else {
        # Run associations against home-made phenotypes. They are assumed to 
        # be quantitative.
        output <- assoc_fast(args[1], args[2])
    }
    write.table(output, args[3], col.names=F, row.names=F, quote=F)
}

if (sys.nframe() == 0) {
    args <- commandArgs(trailingOnly=T)
    print(args)
    main(args)
}
