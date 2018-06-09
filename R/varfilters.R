
# -----------------------------------
# vcfFilter section
# -----------------------------------


#' @title polyIBD vcffilter
#' @description Filtering a variant call file (VCF) using the vcfR package.
#' @param vcffile A variant call file (vcf) path. This VCF will be converted to an object of class \code{vcfR}.
#'
#' @export
#' 

#-----------------------------------------------------
# Read VCF and Go to SNP Matrix of Informative Sites
#------------------------------------------------------
vcffilter <- function(vcffile = NULL, 
                      formatGQ=30, 
                      infoMQ=50,
                      infoQD=25,
                      infoSOR=2,
                      infoAF = 0.05,
                      infoDP = NULL, 
                      prop.loci.missing = 0.05,
                      biallelicSNPs=TRUE){
  
  require(vcfR)
  require(tidyverse)
  
  # -----------------------------------------------------
  # Read and check input
  #------------------------------------------------------
  vcf <- vcfR::read.vcfR(file=vcffile)
  if(biallelicSNPs==TRUE){
    vcf <-vcfR::extract.indels(vcf, return.indels = F) # subset to SNPs
    vcf <- vcf[vcfR::is.biallelic(vcf)] # subset to biallelic
  }else{
    print("polyIBD relies on biallelic SNPs")
  }
  
  # store loci objects on info fields 
  infolist <- c("AF", "DP", "QD", "MQ", "SOR")
  infolist <- lapply(infolist, 
                     function(x){vcfR::extract_info_tidy(vcf, info_fields = x)})
  infodf <- plyr::join_all(infolist, by = "Key", type = "left")
  
  if(typeof(infodf$AF) == "character"){
    infodf$AF <- as.numeric(infodf$AF) # odd default in extract_info -- but character to numeric is fine in R
  }
  
  #--------------------------------------------------------
  # filter loci
  #--------------------------------------------------------
  if(!is.null(infoAF)){
    infodf <- infodf %>% 
      dplyr::filter(AF >= infoAF & AF <= 1-infoAF)
  } 
  if(!is.null(infoDP)){
    DPpercentile <- quantile(infodf$DP, c(infoDP, 1-infoDP))
    infodf <- infodf %>% 
      dplyr::filter(DP >= DPpercentile[1] & DP <= DPpercentile[2])
  }  
  if(!is.null(infoMQ)){
    infodf <- infodf %>% 
      dplyr::filter(MQ >= infoMQ)
  }  
  if(!is.null(infoQD)){
    infodf <- infodf %>% 
      dplyr::filter(QD >= infoQD)
  }  
  if(!is.null(infoSOR)){
    infodf <- infodf %>% 
      dplyr::filter(SOR <= infoSOR)
  }  
  passedloci <- infodf$Key
  
  
  #--------------------------------------------------------
  # filter sample-level GQ
  #--------------------------------------------------------
  # store format and filter fields
  vcfsample <- vcfR::extract.gt(vcf, element =  "GQ", as.numeric = T)
  if(!is.null(formatGQ)){
    vcfsample[vcfsample < formatGQ] <- NA
  }
  vcfsample <- cbind.data.frame(vcf@gt[,1], vcfsample) # need to add format column so matrix vcf@gt will match
  
  #--------------------------------------------------------
  # Subset by loci, GQ
  #--------------------------------------------------------
  vcf@gt[is.na(vcfsample)] <- NA
  
  locimissingness <- apply(vcf@gt, 1, function(x){sum(is.na(x))})
  locimissingness <- which(locimissingness == (ncol(vcf@gt)-1)) # in case a whole loci excluded by GQ
  passedloci <- !passedloci %in% locimissingness
  
  # filter bad loci
  vcf@gt <- vcf@gt[passedloci,]

  #--------------------------------------------------------
  # Drop samples with prop of loci missing
  #--------------------------------------------------------
  sample.prop.loci.missing <- colSums(is.na(vcf@gt))/nrow(vcf@gt)
  vcf@gt <- vcf@gt[, c(prop.loci.missing > sample.prop.loci.missing)]
  
  fix <- as.matrix(vcfR::getFIX(vcf, getINFO = T)[passedloci,])
  gt <- as.matrix(vcf@gt)
  meta <- append(vcf@meta, "##Additional Filters provided by polyIBD filter tools")
  meta <- append(vcf@meta, paste("Some samples may have been filtered by polyIBD filter tools. The new sample count is:", ncol(gt)-1))
  
  # Setting class based off of vcfR documentation https://github.com/knausb/vcfR/blob/master/R/AllClass.R
  newvcfR <- new("vcfR", meta = meta, fix = fix, gt = gt)
  
}


#' @title genautocorr calculates genetic autocorrelation for later linkage disequilibrium filtering
#' @description From an object of class \code{vcfR}, calculate the genetic autocorrelation as the estimate of linkage disequilibrium by genetic distance.
#' @param vcffile A variant call file (vcf) path. This VCF will be converted to an object of class \code{vcfR}.
#'
#' @export
#' 

genautocorr <- function(vcffile = NULL, vcfR = NULL){
  
  require(vcfR)
  require(tidyverse)
  
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
  
  vcf <-vcfR::extract.indels(vcf, return.indels = F) # subset to SNPs
  vcf <- vcf[vcfR::is.biallelic(vcf)] # subset to biallelic
  
  #--------------------------------------------------------
  # Calculate pairwise correlation matrix between observations. Input is a matrix with n rows corresponding to n multivariate measurements, output is a n-by-n correlation matrix. NA values are imputed as the mean.
  #--------------------------------------------------------
  corMat <- function(m) {
    tol <- 1e-9 # tolerance for denominator if AF become the same
    n <- nrow(m)
    c <- matrix(NA,n,n)
    for (i in 1:n) {
      x1 <- unlist(m[i,])
      mu1 <- mean(x1,na.rm=TRUE)
      for (j in i:n) {
        x2 <- unlist(m[j,])
        mu2 <- mean(x2,na.rm=TRUE)
        c[i,j] <- c[j,i] <- sum((x1-mu1)*(x2-mu2),na.rm=TRUE)/sqrt( sum((x1-mu1)^2,na.rm=TRUE)*sum((x2-mu2)^2,na.rm=TRUE) + tol)
      }
    }
    return(c)
  }
  
  # extract the within sample allele frequencies
  vcfAD <- vcfR::extract.gt(vcf, element = "AD")
  vcfAD <- vcfR::AD_frequency(vcfAD)
  vcfAD <- apply(vcfAD, 2, function(x){as.numeric(x)})
  
  # extract distances via positions from vcfR object
  CHROM <- vcfR::getCHROM(vcf)
  POS <- vcfR::getPOS(vcf)
  vcfpos <- cbind.data.frame(CHROM, POS)
  vcfdf <- cbind.data.frame(vcfpos, vcfAD)
  
  vcflist <- split(vcfdf, vcfdf$CHROM)
  
  
  
  cormatgendistwrapper <- function(vcfdf_fromlist){
    
    # get correlation matrix. NA values are imputed as the mean
    df1 <- vcfdf_fromlist[, !colnames(vcfdf_fromlist) %in% c("CHROM", "POS")]
    c <- corMat(df1)
    # get distance between SNPs. This can be extracted from the row names of the vcf
    df2 <- vcfdf_fromlist[, colnames(vcfdf_fromlist) %in% c("POS")]
    gendist <- as.matrix(dist(df2))
    
    ret <- list(vcfAF = vcfdf_fromlist, corMat=c, gendist=gendist)
    return(ret)
  }
  
  
  retlist <- lapply(vcflist, cormatgendistwrapper)
  
  return(retlist)
  
}



#' @title polyIBD filter for linkage disequilibrium 
#' @description Filtering an object of class \code{vcfR} for linkage disequilibrium via genetic autocorrelation.
#' @param vcffile A variant call file (vcf) path. This VCF will be converted to an object of class \code{vcfR}.
#'
#'
#' @export
#' 

genautocorr_filter <- function(vcffile = NULL, vcfR = NULL, genautocorrresult=NULL, threshDist=1e3){
  
  require(vcfR)
  require(tidyverse)
  
  # -----------------------------------------------------
  # Read and check input
  #------------------------------------------------------
  if(is.null(genautocorrresult)){
    stop("Must specify a linkage disequilibrium threshold (see help and tutorial).")
  }
  if(is.null(genautocorrresult)){
    stop("Must specify a genetic autocorrelation results object using the genautocorr function")
  }
  
  if(is.null(vcffile)){
    if(class(vcfR) != "vcfR"){
      stop("vcfR object must be of class vcfR")
    }
    vcf <- vcfR    
  } else{
    vcf <- vcfR::read.vcfR(file=vcffile, verbose=T) # read vcf
  }
  
  vcf <-vcfR::extract.indels(vcf, return.indels = F) # subset to SNPs
  vcf <- vcf[vcfR::is.biallelic(vcf)] # subset to biallelic
  
  vcfdf <- cbind.data.frame(vcf@fix, vcf@gt)
  vcflist <- split(vcfdf, f=factor(vcfdf$CHROM))
  if(length(vcflist) != length(genautocorrresult)){
    stop("The number of chromosomes in the vcfR objec and the results from the genetic autocorrelation analysis differ.")
  }
  # -----------------------------------------------------
  # Filter based on distance
  #------------------------------------------------------
  
  filter_autocorr <- function(vcflist, genautocorrresult){
    gendist <- as.matrix(genautocorrresult$gendist)
    diag(gendist) <- Inf	# block self-comparison
    while (any(gendist<threshDist)) {
      w <- which(gendist<threshDist, arr.ind=TRUE)
      vcflist <- vcflist[-w[1,1],]
      gendist <- gendist[-w[1,1],-w[1,1]]
      
    }
    return(vcflist)
  }
  
  updatedvcflist <- mapply(filter_autocorr, vcflist, genautocorrresult, SIMPLIFY = F)
  updatedvcfdf <- do.call(rbind, updatedvcflist)
  
  # -----------------------------------------------------
  # Return to vcfR object
  #------------------------------------------------------
  fix <- as.matrix(updatedvcfdf[,1:8])
  gt <- as.matrix(updatedvcfdf[,9:ncol(updatedvcfdf)])
  meta <- append(vcf@meta, paste("##Filtered for genetic autocorrelation by polyIBD filter tools with a threshold distance of", threshDist))
  
  # Setting class based off of vcfR documentation https://github.com/knausb/vcfR/blob/master/R/AllClass.R
  newvcfR <- new("vcfR", meta = meta, fix = fix, gt = gt)
  
  return(newvcfR)
  
}