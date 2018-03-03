vcfWindow <-  function(fileName, contig, param, ploidy){
          #read in vcf
          vcf <- readGT(TabixFile(fileName), nucleotide = TRUE, param = param)
          #convert missing to Ns
          vcf[vcf == "."] <- "N/N"
          #make dataframe for wrangling
          vcfDF <- as.data.frame(vcf)
          #separate alleles from either "A/A" or "A|A" to "A" "A"
          windowDF <- lapply(1:ncol(vcfDF), function(x){
            separate(vcfDF, colnames(vcfDF[x]), into = paste(colnames(vcfDF[x]), c(1:ploidy), sep = "/"))[x:(x + (ploidy - 1))]
          }) %>% bind_cols()
          #identify indels as it would ruin the alignement
          indelMat <- lapply(1:ncol(windowDF), function(x){
            nchar(t(windowDF[[x]])) > 1
          })
          #reduce logical vector so that any FALSE == FALSE
          notIndel <- Reduce("|", indelMat)
          #remove from matix and transpose so samples are rows
          alleleMat <- t(windowDF[!notIndel,])
          #make sequence matrix DNAbin from ape package
          as.DNAbin(alleleMat)
}
