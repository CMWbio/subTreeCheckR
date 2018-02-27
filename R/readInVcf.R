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



read.vcf()

winSize <- 100000
nCores <- 5
fileName <- "~/Desktop/Tree-TipR/Plutella_SNPsOnly.vcf.gz"
minSites <- 2000
ploidy <- 2
stat <- c("dist", "pi")
pops <- data.frame(sampleNames = sampleNames, pop = c("PxC", "PxC", "PaC", rep("PxH", 7), "PaG",
                                                      "PxH", "PaR", "PxS", "PaS","PxS", "PaS","PxS",
                                                       "PaS", "PaS", "PxS", "PxG", "PaG", "PaG",
                                                      "PxG", "PaG", "PaC", "PaC", "PaC"))

VCFheader <- scanVcfHeader(fileName)

contigMD <- as.data.frame(VCFheader@header$contig)
contigs <- rownames(contigMD)


prog <- c()
start.time <- Sys.time()
data <- mclapply(contigs[1:10], mc.cores = nCores, function(con){
  prog[1] <<- which(con == contigs)
  length <- as.integer(filter(contigMD, rownames(contigMD) == con)$length)
  if(length >= winSize){
    nWindows <- floor(length / winSize)

    test <- mclapply(seq(1, nWindows), mc.cores = nCores, function(winN){

      pos <- winN * winSize + 1

      p <- ScanVcfParam(which = GRanges(seqnames = con, ranges = IRanges(start = pos - winSize, end = pos)))

      nSites <- length(scanVcf(TabixFile(fileName), param = p)[[1]]$rowRanges)

      if(nSites >= minSites){
        #read in vcf
        vcf <- readGT(TabixFile(fileName), nucleotide = TRUE, param = p)
        #convert missing to Ns
        vcf[vcf == "."] <- "N/N"
        #make dataframe for wrangling
        t <- as.data.frame(vcf)
        #separate alleles from either "A/A" or "A|A" to "A" "A"
        l <- lapply(1:ncol(t), function(x){
          separate(t, colnames(t[x]), into = paste(colnames(t[x]), c(1:ploidy), sep = "/"))[x:(x + (ploidy - 1))]
        }) %>% bind_cols()
        #identify indels as it would ruin the alignement
        p <- lapply(1:ncol(l), function(x){
          nchar(t(l[[x]])) > 1
        })
        #reduce logical vector so that any FALSE == FALSE
        k <- Reduce("|", p)
        #remove from matix and transpose so samples are rows
        s <- t(l[!k,])
        if("pi" %in% stat){
        pi <- lapply(split(pops, pops$pop), function(x){
          #make population matrix
           popGeno <- s[rownames(s) %in% as.vector(outer(as.character(x$sampleNames), 1:ploidy, paste, sep = "/")),]
           # get DNAbin
           dna <- as.DNAbin(popGeno)
           #determine pi for pop
           pi <- data.frame(pi = nuc.div(dna), pop = x$pop[1])
        }) %>% bind_rows() # bind rows together into 1 df
        }
        #convert matrix to DNAbin object from ape
        dna <- as.DNAbin(s)
        #get genetic distance
        if("dist" %in% stat) dist <- dist.dna(dna)

        # to get average distances pull pops out of matrix using vector of names for rows and cols
        #take average and input into data_frame




        # test genind
        # vcf[vcf == "."] <- "N/N"
        # names <- rownames(vcf)
        # names <- gsub("_.*", "", names)
        # names <- gsub("\\.", "_", names)
        # rownames(vcf) <- names
        # genInd <- df2genind(t(vcf), sep = "/", ploidy = 2)
        # dist <- nei.dist(genInd)
        }
      else {
        dist <- NA
      }

      dist
    })
  }
  else{
    test <- "contig too small"
  }
  test
})


end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken





vcf <- vcf[,1:50]
ncol(vcf)


t <- as.data.frame(vcf)

l <- lapply(1:ncol(t), function(x){
  separate(t, colnames(t[x]), into = paste(colnames(t[x]), c(1:2), sep = "/"))[x:(x1)]
}) %>% bind_cols()


test <- dist(t(l), method = "manhattan")


p <- lapply(1:ncol(l), function(x){
  nchar(t(l[[x]])) > 1
})

k <- Reduce("|", p)

s <- t(l[!k,])

test <- as.matrix("Y", "Y")


sequenceNames <- rownames(s)
sampleNames <- gsub("/.*", "", sequenceNames) %>% unique()

dna <- as.DNAbin(s)


dist <- dist.dna(dna)
dist <- as.matrix(dist)
test <- rowMeans2(x = dist, cols = which(sequenceNames == paste(sampleNames[1], 1:2, sep = "/")))

rowMeans2(dist[,1:2])





plot(hclust(dist.dna(dna)))



a <- paste(s[1,], collapse = "") %>% set_names(sampleNames[1])
b <- paste(s[2,], collapse = "") %>% set_names(sampleNames[2])
c <- paste(s[3,], collapse = "") %>% set_names(sampleNames[3])
d <- paste(s[4,], collapse = "") %>% set_names(sampleNames[4])

a <- as.DNAbin(a)

list1 <- list(a,b)
list2 <- list(c,d)
mapply(adist,a,b)

foo <- function(x,y){
  mapply(adist, x, y) / nchar(y))
}
