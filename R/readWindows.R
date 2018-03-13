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
library(d3heatmap)
library(poppr)
library(pegas)
library(phangorn)
library(pbmcapply)

#
winSize <- 50000
nCores <- 5
fileName <- "~/Desktop/Tree-TipR/Plutella_SNPsOnly.vcf.gz"
minSites <- 100
ploidy <- 2
stat <- c("dxy", "pi", "da")
#
# sequenceNames <- rownames(s)
# sampleNames <- gsub("/.*", "", sequenceNames) %>% unique()
# pops <- data.frame(sampleNames = sampleNames, pop = c("PxC", "PxC", "PaC", rep("PxH", 7), "PaG",
#                                                       "PxH", "PaR", "PxS", "PaS","PxS", "PaS","PxS",
#                                                       "PaS", "PaS", "PxS", "PxG", "PaG", "PaG",
#                                                       "PxG", "PaG", "PaC", "PaC", "PaC"))
#
VCFheader <- scanVcfHeader(fileName)

contigMD <- as.data.frame(VCFheader@header$contig)
contigs <- rownames(contigMD)

#
# #alastairu
# pl <- c("KLS0337",
#         "KLS0333","KLS0219","KLS0348","WAM174520",
#         "DARK_TENIUS","KLS0119","MT182","R36639",
#         "MT174","X01049","KLS0114","A00998","X01014",
#         "X01044","X01010","MW04719","KLS0118")
#
# ## Sample/Population dataframe
# pops <- data.frame(sampleNames = pl, pop = c(rep("Afuscus",2), rep("Alaevis",2), rep("Atenius",2), rep("Hcurtus",3),
#                                              "Hcogg", rep("Hcyano",3), rep("Hmelano",2), "Hparvi", rep("Hviper",2)))
#
# contigs <- c("scaffold1|size703937","scaffold2|size540562","scaffold3|size527680",
#              "scaffold4|size481093","scaffold5|size399909","scaffold6|size324431",
#              "scaffold7|size322635","scaffold8|size308202","scaffold9|size294065",
#              "scaffold10|size290244","scaffold11|size286520","scaffold12|size270235")
#


# read in vcf header
vcfHeader <- scanVcfHeader(fileName)

# get contig metadata
contigMD <- as.data.frame(vcfHeader@header$contig)


vcfWindows <- function(fileName, contigs, contigMD, winSize, nCores, ploidy){

  lapply(contigs, function(con){
    length <- as.integer(filter(contigMD, rownames(contigMD) == con)$length)
    if(length >= winSize){
      nWindows <- floor(length / winSize)

     windows <- pbmclapply(seq(1, nWindows), mc.cores = nCores, function(winN){

        pos <- winN * winSize + 1
        start <- pos - winSize
        end <- pos

        p <- ScanVcfParam(which = GRanges(seqnames = con, ranges = IRanges(start = start, end = end)))

        nSites <- tryCatch(length(scanVcf(TabixFile(fileName), param = p)[[1]]$rowRanges),  error=function(e) 0)

        if(nSites >= minSites){
          #read in vcf
          dna <- vcfWindow(fileName = fileName, contig = con, start, end, ploidy = ploidy, haplodize = FALSE, )
          #make sequence matrix DNAbin from ape package
          dna <- list(as.DNAbin(dna))
          names(dna) <- paste0(con, ":", end, "..", start)
        }
        else{
          dna <- list(NA)
          names(dna) <- paste0(con, ":", end, "..", start)
        }
        dna
      }) %>% unlist(recursive = FALSE)
     windows
    }
  }) %>% unlist(recursive = FALSE)

}
