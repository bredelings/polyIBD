#' @title polyIBD vcf2SNPmatrix
#' @description This function wraps the \code{vcfR} package to convert a vcffile or \code{vcfR} object to a SNP matrix suitable for polyIBD. 
#' @param vcffile is the file path that for the VCF. This VCF will be converted to an object of class \code{vcfR}.
#' @param vcfR is the file path that for the VCF. This VCF will be converted to an object of class \code{vcfR}.
#'
#' @export


# TODO:
#     add better error handling
#     Change these functions to one call as an S4 object with slots for SNPmatrix, GT, PopAf, and d
#     better class definitions
#     snpmatrixlist needs to be in an apply loop not for loop 

vcf2infSNPmatrix <- function(vcffile, vcfR) {

  # -----------------------------------------------------
  # Read and check input
  #------------------------------------------------------
  if(is.null(vcffile)){
    if(class(vcfR) != "vcfR"){
      stop("vcfR object must be of class vcfR")
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
  # Determine the number of samples and thier combinations 
  #------------------------------------------------------
  smpls <- factor(colnames(vcfR::extract.gt(vcf)))
  smpls <- t(combn(smpls, 2))
  
  # -----------------------------------------------------
  # Attach Positions to the VCF GT Informative Loci
  #------------------------------------------------------
  
  snpmatrixlist <- lapply(1:nrow(smpls), function(x){matrix(NA, ncol=4, nrow=nrow(snpmatrix))})
    CHROM <- vcfR::getCHROM(vcf)
    POS <- vcfR::getPOS(vcf)
  for(i in 1:nrow(smpls)){
    snpmatrixsave <- snpmatrix[, colnames(snpmatrix) %in% smpls[i,]]
    snpmatrixsave <- cbind.data.frame(CHROM, POS, snpmatrixsave)
    class(snpmatrixsave) <- append(class(snpmatrixsave),"vcfsnpmatrix")
    snpmatrixlist[[i]] <- snpmatrixsave
    rm(snpmatrixsave) # to prevent later clashes
  }
  
  
}


