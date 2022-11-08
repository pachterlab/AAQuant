#!/nfs/fs1/bioinfo/apps-x86_64/R/3.6.1/bin/Rscript
#.libPaths(c("/nfs/odinn/users/kristjanh/lib/Rpackages", .libPaths()))
print(.libPaths())
library(readr)
library(dplyr)
library(tidyr)
library(harmonicmeanp)

aggregate_pvals <- function(pvals_path, weights_path) {
    weights <- read_delim(weights_path, '\t', col_names=c('kmer', 'w'))
    files <- list.files(path=pvals_path, full.names=TRUE, recursive=FALSE)
    output <- data.frame(pheno=unname(sapply(files, basename)), p.hmp=rep(1, length(files)))
    for (file in files) {
        print(file)
        pvals <- read_delim(file, ' ', col_names=c('kmer', 'p'))
        if (all(dim(pvals) == 0)) {
            print("No p values for file")
        } else {
            pvals <- pvals %>% inner_join(weights, by='kmer')
            pval <- tryCatch(
                { p.hmp(pvals$p, w=pvals$w, L=length(pvals$p)) },
                error=function(cond) {
                    return(min(pvals$p))
                },
                warning=function(cond) {
                    return(min(pvals$p))
                }
            )
            output[output$pheno == basename(file),]$p.hmp <- pval
        }
    }  
    return(output)
}      
       
if (sys.nframe() == 0) {
    args <- commandArgs(trailingOnly=T)
    output <- aggregate_pvals(args[1], args[2])
    write.table(output, args[3], sep='\t', row.names=F, col.names=F, quote=F)
}
