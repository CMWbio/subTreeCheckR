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


winSize <- 100000
nCores <- 7
fileName <- "~/Desktop/Tree-TipR/Plutella_SNPsOnly.vcf.gz"
minSites <- 2000
ploidy <- 2
stat <- c("dxy", "pi", "da")

sequenceNames <- rownames(s)
sampleNames <- gsub("/.*", "", sequenceNames) %>% unique()
pops <- data.frame(sampleNames = sampleNames, pop = c("PxC", "PxC", "PaC", rep("PxH", 7), "PaG",
                                                      "PxH", "PaR", "PxS", "PaS","PxS", "PaS","PxS",
                                                      "PaS", "PaS", "PxS", "PxG", "PaG", "PaG",
                                                      "PxG", "PaG", "PaC", "PaC", "PaC"))

VCFheader <- scanVcfHeader(fileName)

contigMD <- as.data.frame(VCFheader@header$contig)

contigsMD <- filter(contigMD, rownames(contigMD) == "KB207950.1") %>% set_rownames("KB207950.1")
contigs <- rownames(contigsMD)

prog <- c()
start.time <- Sys.time()
data <- mclapply(contigs, mc.cores = nCores, function(con){
  length <- as.integer(filter(contigMD, rownames(contigMD) == con)$length)
  if(length >= winSize){
    nWindows <- floor(length / winSize)

    test <- lapply(seq(1, nWindows), function(winN){

      pos <- winN * winSize + 1
      start <- pos - winSize
      end <- pos

      p <- ScanVcfParam(which = GRanges(seqnames = con, ranges = IRanges(start = start, end = end)))

      nSites <- tryCatch(length(scanVcf(TabixFile(fileName), param = p)[[1]]$rowRanges),  error=function(e) 0)

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
        #make sequence matrix DNAbin from ape package
        dna <- as.DNAbin(s)
        #get raw distances using ape::dist.dna as a matrix to calculate dxy, pi and da from
        dist <- dist.dna(dna, as.matrix = TRUE, model = "raw", pairwise.deletion = TRUE)

        #using formula from http://mycor.nancy.inra.fr/egglib/releases/3.0.0a/stats.pdf

        #make pop list
        popList <-  split(pops, pops$pop)

        #calculate dxy for populations from Nei 1987
        dxy <- lapply(popList, function(x){
          lapply(popList, function(y){
            if(!setequal(x$sampleNames, y$sampleNames)){
              #population distance matrix for pairwise pop
              popD <- dist[as.vector(outer(as.character(x$sampleNames), 1:ploidy, paste, sep = "/")), as.vector(outer(as.character(y$sampleNames), 1:ploidy, paste, sep = "/"))]
              #make a tibble with the average number of pairwise differences
              dxy <- data_frame(mean(popD))
              #name col
              colnames(dxy) <-  paste0(x$pop[1], "v" , y$pop[1], "_dxy")
              dxy
            }
          }) %>% bind_cols() #put all onto one row
        })  %>% bind_cols() #bind rows together

        #calculate nucleotide diversity from Nei 1987
        if("pi" %in% stat){
          pi <- lapply(popList, function(x){
            #make population matrix
            popD <- dist[as.vector(outer(as.character(x$sampleNames), 1:ploidy, paste, sep = "/")), as.vector(outer(as.character(x$sampleNames), 1:ploidy, paste, sep = "/"))]
            n <- ncol(popD)
            # get nucleotide diversity
            pi <- sum(popD)/(n*(n-1)/2)
            #determine pi for pop
            piDF <- data_frame(pi)
            # set colnames
            colnames(piDF) <- paste(x$pop[1],"pi", sep = "_")
            piDF
          }) %>% bind_cols() #put all onto one row
        }

        #calculate da from from Nei 1987
        if("da" %in% stat){

         da <- lapply(1:length(dxy), function(x){
           #get dxy pairwise comparison name from colnames of dxy
            dxyName <- colnames(dxy)[1]
            #remove the dxy from the end
            compName <- gsub("_dxy", "", dxyName)
            #get sample names
            samples <- unlist(strsplit(compName, "v"))

            #get pi for population x
            xPi <- pi[[paste0(dxyNameShort[1], "_pi")]]
            #pi for population y
            yPi <- pi[[paste0(dxyNameShort[2], "_pi")]]

            #carry out da calculation from Nei 1987
            da <- dxy[[dxyName]] - ((xPi + yPi)/2)
            #make tibble
            da <- data_frame(da)
            #colnames
            colnames(da) <- paste0(compName, "_da")
            da
          })

        }
      #bind all columns together
      div <- bind_cols(scaffold = con, start = start, end = end, midpoint = (start + end) /2, nSites = nSites, pi, dxy, da)
      }
      else {
      div <- data_frame(scaffold = con, start = start, end = end, midpoint = (start + end) /2, nSites = nSites)
      }
    }) %>% bind_rows() #bind all windows on contig together
  }
  else{
    warning("contig too small")
  }
  test
}) %>% bind_rows() #bind all contigs together


end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken

