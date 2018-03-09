#' Read window from VCF
#'
#' @description Reads in a specified window from a tabix indexed Variant Call Format
#'
#' @details Uses the package \code{VariantAnnotation} to read in a locus (window) from the target VCF,
#'
#'
#' @param fileName A \code{character} vector of length one containing the full path name for a Tabix indexed VCF
#' @param contig A single element \code{character} vector containing the contig name in the VCF
#' @param start object for the target locus
#' @param end end
#' @param ploidy \code{Integer} or \code{Numeric}. The ploidy of the VCF, as called by Variant Caller
#'
#'
#' @return An object of class \code{DNAbin} from the package \code{ape}
#'
#' @importFrom tidyr separate
#'
#' @examples
#'
#' @export
#' @rdname vcfWindow

vcfWindow <-  function(fileName, contig, start, end, ploidy, param, header){

  if(missing(param)) p <- ScanVcfParam(which = GRanges(seqnames = contig, ranges = IRanges(start = start, end = end)))
  else p <- param
  #read in vcf
  vcf <- readGT(TabixFile(fileName), nucleotide = TRUE, param = p)
  #set colnames
  colnames(vcf) <- header@samples

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



}
