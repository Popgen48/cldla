#!/usr/bin/Rscript
library("snpReady")
library("AGHmatrix")

args <- commandArgs(TRUE)


calc_uar_matrix <- function(lgen, output) {
    geno <- read.table(lgen, header = FALSE, na.strings = "NA")
    geno.ready <- raw.data(data = as.matrix(geno), frame = "long", base = TRUE, sweep.sample = 0.5, call.rate = 0.90, maf = 0.025, imput = FALSE)
    M <- geno.ready$M.clean
    geno.ready$report
    Gm <- G.matrix(M = M, method = "UAR", format = "wide")
    write.table(Gm, output, quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)
}

calc_uar_matrix(args[1], args[2])
