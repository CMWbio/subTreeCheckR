gff <- "../Pea Project/References/ref_DBM_FJ_V1.1_scaffolds.gff3"
fileName <- "~/Desktop/Tree-TipR/Plutella_SNPsOnly.vcf.gz"

getGenes <- function(fileName, gff, feature = "gene", featureCoverage = 0.1, nCores = 5, fastaOut = "~/Desktop/Tree-TipR/fasta_test/"){

  allGR <- import.gff(gff)


  if(feature == "gene"){

    genes <- allGR[allGR$type == "gene"]
    featureList <- split(genes, genes$Name)

  }

  if(feature == "exon"){

    exons <- allGR[allGR$type == "exon"]
    featureList <- split(exons, unlist(exons$Parent))

    featureList <- mclapply(1:length(featureList), mc.cores = nCores, function(x){

      gr <- featureList[[x]]

      logi <- allGR$ID == gr$Parent[[1]]
      logi[is.na(logi)] <- FALSE

      gr$Name <- allGR[logi]$Name

      gr
    })

  }

  vcfHeader <- scanVcfHeader(fileName)

  genesList <- pbmclapply(1:length(featureList), mc.cores = nCores, function(x){

    p <- ScanVcfParam(which = featureList[[x]])

    minSites <- round(sum(featureList[[1]]@ranges@width) * featureCoverage)

    nSites <- tryCatch(length(scanVcf(TabixFile(fileName), param = p)[[1]]$rowRanges),  error=function(e) 0)
    if(nSites >= minSites){
      dat <- vcfWindow(fileName = fileName, ploidy = 2, param = p, header = vcfHeader)
      dna <- as.DNAbin(dat)

      dnaCat <- lapply(dna, function(x) as.character(x[1:length(x)]))

      if(feature == "gene") {
        geneName <- featureList[[x]]$Name[1]}else{
          geneName <- paste0(featureList[[x]]$Name[1], "_exons")}

      fastaName <- paste0(fastaOut, "/", geneName)
      cat(file=fastaName, paste(paste0(">",names(dnaCat)), sapply(dnaCat, paste, collapse=""), sep="\n"), sep="\n")
      dna
    }
  })
  genesList
}
