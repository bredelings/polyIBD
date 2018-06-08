#' @title polyIBD vcf2polyIBDinput
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
#     better management of memory with p

vcf2polyIBDinput <- function(vcffile=NULL, vcfR=NULL) {

  # -----------------------------------------------------
  # Read and check input
  #------------------------------------------------------
  if(is.null(vcffile)){
    if(class(vcfR) != "vcfR"){
      stop("vcfR object must be of class vcfR")
    }
    vcf <- vcfR
  } else if (!is.null(vcffile)){
    vcf <- vcfR::read.vcfR(file=vcffile, verbose=F) # read vcf
  } else {
    stop("Must specify an input")
  }
  
  vcf <-vcfR::extract.indels(vcf, return.indels = F) # subset to SNPs
  vcf <- vcf[vcfR::is.biallelic(vcf)] # subset to biallelic
  print(vcf)
  # -----------------------------------------------------
  # determine ploidy to determine genotype numeric placeholder  
  #------------------------------------------------------
  ploidy <- stringr::str_extract(vcf@meta[grepl("ploidy", vcf@meta)], "ploidy=[0-9]")
  ploidy <- stringr::str_split(t(ploidy), "=", simplify = T)
  ploidy <- as.numeric(ploidy[1,2])

  if(ploidy == 2){
    
    snpmatrix <- vcfR::extract.gt(vcf, element='GT', as.numeric=F) # numeric as T doesn't parse 0/1 correctly
    snpmatrix[snpmatrix == "0/0"] <- 0
    snpmatrix[snpmatrix == "0/1"] <- 1
    snpmatrix[snpmatrix == "1/1"] <- 2
    snpmatrix <- apply(snpmatrix, 2, function(x){as.numeric(x)}) # need to convert from char (--dependent on case of "/") to numeric
  } else if(ploidy == 1){
    snpmatrix <- vcfR::extract.gt(vcf, element='GT', as.numeric=T)
    snpmatrix[snpmatrix == 0] <- 0
    snpmatrix[snpmatrix == 1] <- 2
  } else {
    stop("You have a ploidy that is less than 1 or greater than 3, which cannot be accomodated by polyIBD")
  }
  # -----------------------------------------------------
  # Determine the number of samples and thier combinations 
  #------------------------------------------------------
  smpls <- factor(colnames(vcfR::extract.gt(vcf)))
  smpls <- t(combn(smpls, 2))
  
  CHROM <- vcfR::getCHROM(vcf)
  POS <- vcfR::getPOS(vcf)
  
  # -----------------------------------------------------
  # Calculate pop AF, p
  #------------------------------------------------------
  L <- length(POS)
  p <- rep(NA, length(POS))
  p <- apply(snpmatrix, 1,
             function(x){(2*(length(which(x==0))) + length(which(x==1)))/(2*length(x))}) # since we know homozygous ref is 0, so this counts as 2 As, and then we count hets and then divide by 2*possible alleles. Doing this roundabout way because we aren't in HWE
  
  
  # -----------------------------------------------------
  # Attach Positions to the VCF GT Informative Loci
  #------------------------------------------------------
  
  retlist <- lapply(1:nrow(smpls), function(x){
    list(samples = rep(NA, 1),
         snpmatrix = matrix(NA, ncol=4, nrow=nrow(snpmatrix)),
         p=rep(NA, nrow(snpmatrix)))
    })
    
  for(i in 1:nrow(smpls)){
    snpmatrixsave <- snpmatrix[, colnames(snpmatrix) %in% smpls[i,]]
    snpmatrixsave <- cbind.data.frame(CHROM, POS, snpmatrixsave)
    retlist[[i]]$snpmatrix <- snpmatrixsave
    retlist[[i]]$p <- p
    retlist[[i]]$samples <- paste(smpls[i,], collapse = "||")
    
  }
  
  # -----------------------------------------------------
  # return
  #-----------------------------------------------------
  
  class(retlist) <- "polyIBDinput"
  return(retlist)
  
}


