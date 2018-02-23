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
#' @return UI data for fastQC shiny.
#'
#' @examples
#' \dontrun{
#'  samples <- rownames(distMList[[1]])[order(c(rownames(distMList[[1]])))]
#'  groups <-  c(rep("I2", 4), "t", "t", rep("I1", 4), "t", "t", rep("O", 4))
#'
#' @export
#' @rdname threePop
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
library(reshape2)
library(phangorn)
library(phytools)

subTrees <- read.tree("../PANAMA.ALLSHARED.w100.rerooted.nwk")

distMList<- subTrees %>% lapply(FUN = function(x){
  matrix <- cophenetic(x)
})

samples <- rownames(distMList[[1]])[order(c(rownames(distMList[[1]])))]
groups <-  c(rep("P3", 4), "t", "t", rep("P2", 4), "t", rep("P1", 4), "t")
metaD <- data_frame(samples = samples, groups = groups)



threePop <- function(subTrees, metaD, cores = 1){

  stopifnot(is.data.frame(metaD))
  if(length(metaD) != 2 | !all(colnames(metaD) == c("samples", "groups"))) stop("Object metaD does not contain 2 columns named 'samples' and 'groups'")
  if(!all(c("P1", "P2", "P3") %in% unique(metaD$groups))) stop("Column 'groups' metaD object does not contain the groups P1, P2 and P3")
  if(!is.numeric(cores) | cores < detectCores()) stop("Value of cores is either non-numeric or greater than total number of cores - 1")
  if(class(subTrees) == "character") subTrees <- read.tree(subTrees)

  distMList<- subTrees %>% lapply(FUN = function(x){
    matrix <- cophenetic(x)
  })

  metaDList <- split(metaD, groups)
  P1 <- metaDList$P1[[1]]
  P2 <- metaDList$P2[[1]]
  P3 <- metaDList$P3[[1]]

  testAll <- mclapply(seq(length(distMList)), mc.cores = cores, function(window){
    win <- distMList[[window]]
    lapply(P1, function(x){
      lapply(P3, function(y){
        nullDistance <-  mean(unlist(lapply(P2, function(z) {win[z, y]})))
        testDistance <- win[x, y]
        data_frame(TSP = testDistance/(nullDistance+testDistance), pairwiseComp = paste(x, y, sep =":"), window = window)
      }) %>% bind_rows()
    }) %>% bind_rows()
  }) %>% bind_rows()
}
