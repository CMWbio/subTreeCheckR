#' Three Population Tree Tip Proportion
#'
#' @description Carries out Three Population Subtree Tip proportions from \code{Ward & Baxter (2018), GBE}
#'
#' @details Currently some plots can take a while to render if the \code{FastqcDataList} passed to
#' \code{fastqcInput} has many elements
#'
#' @param subTrees Trees to carry out 3 populaiton Tree Tip proportions on. Can be either a \code{character} file path or \code{multiphy}
#' object from the package \code{ape}
#' @param metaD Trees to carry out 3 populaiton Tree Tip proportions on. Can be either a \code{character} file path or \code{multiphy}
#' object from the package \code{ape}
#'
#'
#' @return \code{list} of trees in \code{phy} format.
#'
#' @examples
#'
#' @export
#' @rdname buildWindowTrees
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


buildWindowTrees <- function(fileName, contigs, winSize = 100000,
                          subModel = "raw", minSites, ploidy = 2,
                          output = "tree", nCores = 2){

  # checks
  if(!is.character(filename) | length(fileName) > 1) stop("fileName must ba a character vector of length 1")
  stopifnot(file.exists(fileName))
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
  if(missing(minSites)) minSites <- 1
  # if contigs is missing get contigs from contigMD (all contigs)
  if(missing$contigs) contigs <- rownames(contigMD)



  nestedList <- lapply(contigs, function(con){
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

        dna <- vcfWindow(fileName = fileName, contig = con, param = p, ploidy = ploidy)

        div <- list(dist.dna(dna, model = subModel, pairwise.deletion = TRUE))
        names(div) <- paste0(con, ":", (end - start)/2)
        }
        else {
          div <- NA
        }
        div
      }) #bind all windows on contig together
    }
    else{
      warning("contig too small")
    }
    distanceList
  }) %>% unlist()
}
