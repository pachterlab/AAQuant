#
# Rscript --vanilla filter_and_transform_phenotypes.R /nfs/odinn/users/kristjanh/work/kmer_association/dbg/prune_pipeline/assoc/phenotype_lists/quantitative_phenotypes /nfs/odinn/users/kristjanh/work/kmer_association/dbg/prune_pipeline/pn_lists/blood_bamfiles.filtered.pns /nfs/odinn/users/kristjanh/work/kmer_association/dbg/prune_pipeline/assoc/phenotypes_05052020
#

library(readr)
library(dplyr)
library(tidyr)

rank_transform <- function(x){
    return(qnorm((rank(x,ties.method='random')-0.5)/length(x)))
}

main <- function(args) {

    pns_path <- "/nfs/odinn/users/kristjanh/work/kmer_association/dbg/prune_pipeline/pn_lists/blood_bamfiles.filtered.pns"
    pheno_path     <- args[1]
    pns_path       <- args[2]
    output_path    <- args[3]
    threshold_path <- args[4]

    pns <- read_lines(pns_path)

    paths <- read_delim(pheno_path, '\t', col_names=c('pheno', 'path'))

    for (i in 1:dim(paths)[1]) {

        pheno <- paths[i,'pheno']
        path  <- paths[i,'path']

        if (file.exists(file.path(path, 'qtl.val'))) {
            phenos <- read_delim(file.path(path, 'qtl.val'), delim='\t', col_names=c('PN', 'pheno'))
        } else if (file.exists(file.path(args[2], 'work/qtl.raw'))) {
            phenos <- read_delim(file.path(path, 'work/qtl.raw'), delim='\t', col_names=c('PN', 'pheno'))
        } else {
            print('File ' + path + ' not found.')
            next
        }
        out <- phenos[phenos$PN %in% pns,]
        out$pheno <- rank_transform(out$pheno)

        threshold <- qt((1 - (0.05/20000)/2), dim(out)[1] - 2)

        write_delim(out, file.path(output_path, pheno), append=F, col_names=F, quote_escape=F)
        cat(threshold, file=file.path(threshold_path, pheno), append=F)
    }

}

if (sys.nframe() == 0) {
    args <- commandArgs(trailingOnly=T)
    print(args)
    main(args)
}
