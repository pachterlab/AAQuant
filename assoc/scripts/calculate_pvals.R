library(readr)
library(dplyr)
library(tidyr)

calculate_pvals <- function(path) {
    assoc <- read_delim(path, delim='\t', col_names=T)
    assoc$p <- 2*pt(abs(assoc$coef/assoc$stderr), assoc$df, lower.tail=F)
    write_delim(assoc, path, delim='\t', append=F, col_names=T, quote_escape=F)
}

main <- function(args) {
    dirs <- list.files(args[1], full.names=T, recursive=F, include.dirs=T)
    for (d in dirs) {
        files <- list.files(file.path(d, 'assoc'), full.names=T, recursive=F, include.dirs=F)
        for (f in files) {
            calculate_pvals(f)
        }
    }
}

if (sys.nframe() == 0) {
    args <- commandArgs(trailingOnly=T)
    main(args)
}
