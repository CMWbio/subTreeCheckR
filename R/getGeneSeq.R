getGeneSeq <- function(fileName, genes, fastaOut = "~/Desktop/Tree-TipR/fasta_test/",
                       haploidize = TRUE, ploidy = 2, featureCoverage = 0.1, feature = "gene", nCores = 1,
                       minExp = 0.5){

  vcfHeader <- scanVcfHeader(fileName)

  genesList <- mclapply(1:length(genes), mc.cores = nCores, function(x){

    p <- ScanVcfParam(which = genes[[x]])

    minSites <- round(sum(genes[[1]]@ranges@width) * featureCoverage)

    nSites <- tryCatch(length(scanVcf(TabixFile(fileName), param = p)[[1]]$rowRanges),  error=function(e) 0)
    if(nSites >= minSites){
      dat <- vcfWindow(fileName = fileName, ploidy = ploidy, param = p, header = vcfHeader, haploidize = haploidize)

      if(all(rowSums(dat == "N")/ncol(dat) < minExp)){
        dna <- as.DNAbin(dat)

        if(!is.na(fastaOut)){
          dnaCat <- lapply(dna, function(x) as.character(x[1:length(x)]))

          if(feature == "gene") {
            geneName <- genes[[x]]$Name[1]}else{
              geneName <- paste0(genes[[x]]$Name[1], "_exons")}

          fastaName <- paste0(fastaOut, "/", geneName)
          cat(file=fastaName, paste(paste0(">",names(dnaCat)), sapply(dnaCat, paste, collapse=""), sep="\n"), sep="\n")

        }

        dna <- list(dna)
        names(dna) <- genes[[x]]$Name[1]
        dna
      }
      else{
        NULL

      }

    }

  }) %>% unlist(recursive = FALSE)

  genesList
}
