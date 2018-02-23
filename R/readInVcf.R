library(readr)
library(ape)
library(magrittr)
library(dplyr)
library(tidyr)
library(scales)
library(parallel)
library(ggplot2)
library(plotly)
library(geiger)
library(phytools)
library(VariantAnnotation)
library(stringr)
library(msa)




p <- ScanVcfParam(which = GRanges(seqnames = "1", ranges = IRanges(start = 1, end = 500000)))
mclapply(1:1000, mc.cores = 5, function(x){
vcf <- readGT(TabixFile("~/Downloads/1.1-1000000.ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes(1).vcf.gz"), nucleotides = TRUE,
              param = p)
NULL

})

ncol(vcf)
vcf <- vcf[1:50,]

t <- as.data.frame(vcf)

l <- lapply(1:ncol(t), function(x){
  separate(t, colnames(t[x]), into = paste(colnames(t[x]), c(1:2), sep = "."))[x:(x+1)]
}) %>% bind_cols()

p <- lapply(1:ncol(l), function(x){
nchar(t(l[[x]])) > 1
})

k <- Reduce("|", p)




s <- t(l[!k,])

sampleNames <- rownames(s)

a <- paste(s[1,], collapse = "") %>% set_names(sampleNames[1])
b <- paste(s[2,], collapse = "") %>% set_names(sampleNames[2])
mapply(adist,a,b)


