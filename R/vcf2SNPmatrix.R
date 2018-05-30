#' @title polyIBD vcf2SNPmatrix
#' @description Wrapper of the vcfR package.  Goal of readVcfsplittoSNPMATRIXlist is to take a vcf, convert it to a matrix, find informative sites. Goal of getLociPopAF is to calculate population allele frequency by loci. 
#' @param vcffile is the file path that for the VCF. This VCF will be converted to an object of class \code{vcfR}.
#'@param vcfR is the file path that for the VCF. This VCF will be converted to an object of class \code{vcfR}.
#'
#' @export


# TODO:
#     add better error handling
#     Change these functions to one call as an S4 object with slots for SNPmatrix, GT, PopAf, and d


vcf2infSNPmatrix <- function(vcffile, vcfR) {

  # -----------------------------------------------------
  # Read and check input
  #------------------------------------------------------
  if(is.null(vcffile)){
    if(class(vcfR) != "vcfR"){
      paste("vcfR object must be of class vcfR")
    }
    vcf <- vcfR
  } else{
    vcf <- vcfR::read.vcfR(file=vcffile, verbose=T) # read vcf
  }
  
  # -----------------------------------------------------
  # determine ploidy to determine genotype numeric placeholder  
  #------------------------------------------------------
  ploidy <- stringr::str_extract(vcf@meta[grepl("ploidy", vcf@meta)], "ploidy=[0-9]")
  ploidy <- stringr::str_split(t(ploidy), "=", simplify = T)
  ploidy <- as.numeric(ploidy[1,2])

  if(ploidy == 2){
    
    snpmatrix <- vcfR::extract.gt(vcf, element='GT', as.numeric=F)
    snpmatrix[snpmatrix == "0/0"] <- 0
    snpmatrix[snpmatrix == "0/1"] <- 1
    snpmatrix[snpmatrix == "1/1"] <- 2
    
  } else if(ploidy == 1){
    snpmatrix <- vcfR::extract.gt(vcf, element='GT', as.numeric=F)
    snpmatrix[snpmatrix == "0"] <- 0
    snpmatrix[snpmatrix == "1"] <- 2
    
  } else {
    print("You have a ploidy that is less than 1 or greater than 3, which cannot be accomodated by polyIBD")
  }
  
  # -----------------------------------------------------
  # Attach Positions to the VCF GT Informative Loci
  #------------------------------------------------------

  REG <- data.frame(datvcf@fix, stringsAsFactors = F)
  REG <- REG[,colnames(REG) %in% c("CHROM", "POS")] #vcf format it is first two columns but to be explicit
  REG$POS <- as.numeric(REG$POS) # taking char to num -- fine in R, maybe not in Cpp
  snpmatrix <- cbind(REG, snpmatrix)
  return(snpmatrix)
}


