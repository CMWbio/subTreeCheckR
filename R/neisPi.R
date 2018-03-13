neisPi <- function(dist, popList){

  lapply(popList, function(x){
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
