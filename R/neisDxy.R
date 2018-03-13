neisDxy <- function(dist, popList) {

  dxy <- lapply(popList, function(x){
  lapply(popList, function(y){
    if(!setequal(x$sampleNames, y$sampleNames)){
      #population distance matrix for pairwise pop
      popD <- dist[as.vector(outer(as.character(x$sampleNames), 1:ploidy, paste, sep = "/")), as.vector(outer(as.character(y$sampleNames), 1:ploidy, paste, sep = "/"))]
      #make a tibble with the average number of pairwise differences
      dxy <- data_frame(mean(popD), sd(popD))
      #name col
      colnames(dxy) <-  paste0(x$pop[1], "v/" , y$pop[1], c("_dxy", "_SDdxy"))
      dxy
    }
  }) %>% bind_cols() #put all onto one row
})  %>% bind_cols() #bind rows together

}
