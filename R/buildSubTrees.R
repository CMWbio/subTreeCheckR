#' Build a tree for each window
#'
#' @description Build local trees for each window along a contig.
#'
#' @details Builds a tree using data from each window across a contig. Defined as a local tree these
#' represent the evolutionary history of the window.
#'
#' @param fileName \code{character}. Full path name for a Tabix indexed VCF file.
#' @param DNAwin \code{list} of multiple \code{DNAbin}. A list of already imported VCF files.
#' @param contigs \code{character}. Default is \code{"all"} Contigs to extract windows from.
#' @param winSize \code{numeric}. Default is \code{100000}. Window size in base pairs.
#' @param subModel \code{character}. The substitution model to fit to the data. Accepts nucleotide models defined in \code{ape::dist.dna}:
#' "raw", "N", "TS", "TV", "JC69", "K80", "F81", "K81", "F84", "BH87", "T92", "TN93", and "GG95."
#' @param minSites \code{numeric}. Default is 10 percent of winSize. The minimum number of sites within a window to consider it.
#' @param ploidy \code{numeric}. The ploidy of the VCF, as called by Variant Caller.
#' @param nCores \code{numeric}. Number of cores to run analysis on.
#' @param write \code{character}. Default is NA, will only return the R object. The directory to write the trees out to.
#'
#'
#' @return \code{list} of trees in \code{phy} format. Newick format trees in a specified directory.
#'
#' @import VariantAnnotation
#' @import ape
#' @importFrom parallel mclapply
#' @importFrom dplyr filter
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#'
#' @examples
#'
#' @export
#' @rdname buildLocalTrees
# library(readr)
# library(ape)
# library(magrittr)
# library(dplyr)
# library(tidyr)
# library(scales)
# library(parallel)
# library(ggplot2)
# library(plotly)
# library(geiger)
# library(phytools)
# library(VariantAnnotation)
# library(stringr)
# library(msa)
# library(d3heatmap)
# library(poppr)
# library(pegas)
# library(phangorn)


buildLocalTrees <- function(fileName, DNAwin = NA, contigs, winSize = 100000,
                            subModel = "raw", minSites, ploidy = 2,
                            nCores = 5, write = NA, haploidize = TRUE){

  # checks
  if(!is.character(fileName) | length(fileName) > 1) stop("fileName must ba a character vector of length 1")
  stopifnot(file.exists(fileName))
  # stopifnot(is.list(DNAwin))
  if(!missing(contigs) && !is.character(contigs)) stop("fileName must ba a character vector of length 1")
  if(!is.numeric(winSize)) stop("winSize should be numeric")
  if(!is.character(subModel) | length(subModel) >1) stop("subModel must be a character vector of length 1")
  if(!subModel %in% c("raw", "N", "TS", "TV", "JC69", "K80", "F81",
                      "K81", "F84", "BH87", "T92", "TN93", "GG95", "logdet", "paralin", "indel", "indelblock"))
    if(!type %in% c("tree", "distance matrix"))
      stopifnot(nCores < detectCores() - 1)



  # read in vcf header
  vcfHeader <- scanVcfHeader(fileName)

  # get contig metadata
  contigMD <- as.data.frame(vcfHeader@header$contig)

  # set minSites to 1
  if(missing(minSites)) minSites <- 0.05 * winSize
  # if contigs is missing get contigs from contigMD (all contigs)
  if(missing(contigs)) contigs <- rownames(contigMD)


  if(any(is.na(DNAwin))){
  nestedList <- pblapply(contigs, function(con){
    length <- as.integer(filter(contigMD, rownames(contigMD) == con)$length)
    if(length >= winSize){
      nWindows <- floor(length / winSize)
        distanceList <- mclapply(seq(1, nWindows), mc.cores = nCores, function(winN){

          pos <- winN * winSize + 1
          start <- pos - winSize
          end <- pos

          p <- ScanVcfParam(which = GRanges(seqnames = con, ranges = IRanges(start = start, end = end)))

          nSites <- tryCatch(length(scanVcf(TabixFile(fileName), param = p)[[1]]$rowRanges),  error=function(e) 0)

          if(nSites >= minSites){

            dna <- vcfWindow(fileName = fileName, contig = con, param = p, ploidy = ploidy, haploidize = haploidize, header = vcfHeader)
            dna <- as.DNAbin(dna)
            dist <- dist.dna(dna, model = "K80", pairwise.deletion = TRUE)
            tree <- upgma(dist)
            tree$edge.length <- abs(tree$edge.length)
            tree <- list(tree)
            names(tree) <- paste0(con, ":", end, "..", start)

          }
          else {
            tree <- list(NA)
          }
          tree
        }) %>% unlist(recursive = FALSE)#bind all windows on contig together
      }
    else{
      distanceList <- list(NA)
    }
    distanceList
  }) %>% unlist(recursive = FALSE)
 }
  else {
    nestedList <- lapply(1:15, function(x){
      dist <- dist.dna(DNAwin[[x]], model = subModel, pairwise.deletion = TRUE)
      tree <- upgma(dist)
      tree$edge.length <- abs(tree$edge.length)
      tree <- list(tree)
      names(tree) <- names(DNAwin)[x]
      tree
    }) %>% unlist(recursive = FALSE)
  }

  nestedList <- nestedList[!is.na(nestedList)]
  class(nestedList) <- "multiPhylo"

  if(!is.na(write)) write.tree(phy = nestedList, tree.names = TRUE, file = write)

  nestedList
}


