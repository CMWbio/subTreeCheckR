#' Three Population Tree Tip Proportion
#'
#' @description Carries out Three Population Subtree Tip proportions from \code{Ward & Baxter (2018), GBE}
#'
#' @details Currently some plots can take a while to render if the \code{FastqcDataList} passed to
#' \code{fastqcInput} has many elements
#'
#' @param localTrees Trees to carry out 3 populaiton Tree Tip proportions on. Can be either a \code{character} file path or \code{multiphy}
#' object from the package \code{ape}
#' @param metaD Trees to carry out 3 populaiton Tree Tip proportions on. Can be either a \code{character} file path or \code{multiphy}
#' object from the package \code{ape}
#' @param nnCores \code{numeric}. Number of nCores to run analysis on.
#'
#'
#' @return UI data for fastQC shiny.
#'
#'
#' @examples
#'
#' @export
#' @rdname threePop
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
# library(reshape2)
# library(phangorn)
# library(phytools)
#
localTrees <- read.tree("~/Desktop/Tree-TipR/PANAMA.ALLSHARED.w100.rerooted.nwk")
#
# localTrees(localTrees[[1]])
#
#
#
#
# distMList<- localTrees %>% lapply(FUN = function(x){
#   matrix <- cophenetic(x)
# })
#
samples <- rownames(distMList[[1]])[order(c(rownames(distMList[[1]])))]
groups <-  c(rep("P3", 4), "t", "t", rep("P2", 4), "t", rep("P1", 4), "t")
 metaD <- data_frame(samples = samples, groups = groups)
# #


threePop <- function(localTrees, genomeTree, metaD, nCores){

  if(missing(genomeTree)) genomeTree <- consensus(localTrees, p = 0.5, check.labels = TRUE)

  stopifnot(is.data.frame(metaD))
  if(length(metaD) != 2 | !all(colnames(metaD) == c("samples", "groups"))) stop("Object metaD does not contain 2 columns named 'samples' and 'groups'")
  if(!all(c("P1", "P2", "P3") %in% unique(metaD$groups))) stop("Column 'groups' metaD object does not contain the groups P1, P2 and P3")
  if(!is.numeric(nCores) | nCores > detectCores()) stop("Value of nCores is either non-numeric or greater than total number of nCores - 1")
  if(class(localTrees) == "character") localTrees <- read.tree(localTrees)

  distMList<- localTrees %>% lapply(FUN = function(x){
    matrix <- cophenetic(x)
  })

  metaDList <- split(metaD, groups)
  P1 <- metaDList$P1[[1]]
  P2 <- metaDList$P2[[1]]
  P3 <- metaDList$P3[[1]]

  testAll <- mclapply(seq(length(localTrees)), mc.cores = 4, function(window){
    win <- distMList[[window]]
    lapply(P1, function(x){
      lapply(P3, function(y){
        nullDistance <-  mean(unlist(lapply(P2, function(z) {win[z, y]})))
        testDistance <- win[x, y]
        data_frame(TSP = testDistance/(nullDistance+testDistance), pairwiseComp = paste(x, y, sep =":"), window = window)
      }) %>% bind_rows()
    }) %>% bind_rows()
  }) %>% bind_rows()


  distMgen <- cophenetic.phylo(genomeTree)

  testGen <- lapply(P1, function(x){
    lapply(P3, function(y){
      nullDistance <-  mean(unlist(lapply(P2, function(z) {distMgen[z, y]})))
      testDistance <- distMgen[x, y]
      data_frame(TSP = testDistance/(nullDistance+testDistance), pairwiseComp = paste(x, y, sep =":"), window = window)
    }) %>% bind_rows()
  }) %>% bind_rows()

  splitWi <- split(testAll, testAll$pairwiseComp)
  splitWg <- split(testGen, testGen$pairwiseComp)

  statW <- mclapply(splitWi, mc.cores = nCores, function(x){

    genComp <- splitWg[[x$pairwiseComp[1]]]

    logit <- function(x){

      -log((1/x) - 1)

    }

    t <- logit(x$TSP) - logit(genComp$TSP)
    u <- mean(t)


    reSam <- lapply(1:500, function(x){
      n <- round(0.01 * length(t))
      sim <- sample_n(data_frame(t), size = n)
      u_sim <- mean(sim$t)
    }) %>% unlist()

    SD <- sd(reSam)
    Z <- u/SD

    data_frame(compName = x$pairwiseComp[1], mean = u, "Jackknife SE" =  SD, Z = Z)

  }) %>% bind_rows()

  statW
}






#
# start.time <- Sys.time()
# t1 <- threePop(localTrees = localTrees, metaD = metaD)
# end <- Sys.time()
# end - start.time
