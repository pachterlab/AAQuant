#!/nfs/fs1/bioinfo/apps-x86_64/R/3.6.1/bin/Rscript
print(.libPaths())
print(R.version.string)
library(readr)
library(dplyr)
library(tidyr)

rank_transform <- function(x){
    return(qnorm((rank(x,ties.method='random')-0.5)/length(x)))
}

assoc_quant <- function(expr_path, pheno_path) {
    df_pheno <- read_delim(pheno_path, '\t', col_names=c('PN', 'gt'))
    df <- read_delim(expr_path, '\t', col_names=c('PN', 'kmer', 'expr'))
    pns <- pull(arrange(distinct(df, PN), PN))
    df_pheno <- filter(df_pheno, PN %in% pns)
    df_pheno$gt <- rank_transform(df_pheno$gt)

    df <- df %>% complete(kmer, nesting(PN, gt), fill=list(expr=0))
    df <- pivot_wider(df, id_cols=PN, names_from=kmer, values_from=expr)
    df <- df %>% inner_join(df_pheno, by='PN')

    kmers <- pull(arrange(distinct(df, kmer), kmer))

    form <- as.formula(paste('gt', '~', paste(kmers, collapse=' + '), '-1', collapse=' '))
    ll <- glm(form, family=gaussian(), data=df)
    cc <- coef(summary(ll))
    if (length(cc[,'Pr(>|t|)']) >= 1) {
        output$p[output$kmer == km] <- cc['expr','Pr(>|t|)']
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
        output <- assoc_quant(args[1], args[2])
    }
    write.table(output, args[3], col.names=F, row.names=F, quote=F)
}

if (sys.nframe() == 0) {
    args <- commandArgs(trailingOnly=T)
    print(args)
    main(args)
}
