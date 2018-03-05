#' Compare trees to clade structure
#'
#' @description Identifies trees discordant with a tree topology
#'
#' @details checks a all trees for a provided clade structure
#'
#'
#' @param trees \code{multiphy} or \code{list} of \code{phy}.
#' Trees to test for clade structure specified in \code{cladeStr}
#' @param cladeStr \code{data.frame}. Clade structure of the tree.
#'
#'
#' @return \code{list} of trees in \code{phy} format containing all discordant trees.
#'
#' @importFrom phytools getDescendants
#'
#' @examples
#'
#'
#' @export
#' @rdname topologyCheck



# library(reeadr)
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
#
# cladeStr <- data_frame(tipLabels = , cladeLabels = c(rep("I1", 8),rep("I2", 8),rep("O", 8)))
#
# trees <- read.tree("../Tree-TipR/PANAMA.ALLSHARED.w100.rerooted.nwk")
#
# subTrees <- read.tree("../PANAMA.ALLSHARED.w100.rerooted.nwk")
#
# distMList<- trees %>% lapply(FUN = function(x){
#   matrix <- cophenetic(x)
# })
#
# samples <- rownames(distMList[[1]])[order(c(rownames(distMList[[1]])))]
# groups2 <-  c(rep("P3", 4), "t", "t", rep("P2", 4), "t", rep("P1", 4), "t")
# groups <-  c(rep("P3", 4), "t", "t", rep("P2", 4), "t", rep("P2", 4), "t")
# cladeStr <- data_frame(samples = samples, cladeLabels = groups, groups = groups2)
#
#
# clad <- test[[1]]

#####

topologyCheck <- function(trees, cladeStr){
  allTrees <- mclapply(seq(length(trees)), mc.cores = 5, function(treeN){
    tree <- trees[[treeN]]
    data <- c()
      for (x in 2:length(cladeStr)) {
      data[x-1] <- split(cladeStr, cladeStr[x]) %>% lapply(function(clad){
      # get target taxa from metaData
      targetTaxa <- clad[[1]]
      #get the number of nodes in the Subtree
      Nnodes <- tree[["Nnode"]]
      #get all tips in tree
      tips <- tree[["tip.label"]]
      #get the MRCA of targetTaxa
      MRCA <- getMRCA(tree, targetTaxa)
      #get the clade descended from the the MRCA
      clade <- getDescendants(tree, MRCA)
      #total number of tips in tree
      tipN <- length(tips)
      #remove internal nodes from the clade
      cladeTips <- clade[clade <= tipN]
      #get names for tips in clade
      cladeTipNames <- tips[c(1:tipN) %in% cladeTips]
      setequal(targetTaxa, cladeTipNames)

    }) %>% unlist() %>% all()
      if(!any(data)) break
      }
    if(!data){
      data.frame(treeType = "Discordant", window = treeN)
    }

      }) %>% bind_rows()
  }









