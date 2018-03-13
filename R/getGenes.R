gff <- "../Pea Project/References/ref_DBM_FJ_V1.1_scaffolds.gff3"
fileName <- "~/Desktop/Tree-TipR/Plutella_SNPsOnly.vcf.gz"

getGenes <- function(fileName, gff, feature = "gene", nCores = 5){

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
  featureList[1:100]


}
