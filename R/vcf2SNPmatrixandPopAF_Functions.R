#' @title polyIBD Transition Probabilities 
#' @description .....
#' @param file
#' @export

#   Goal of readVcfsplittoSNPMATRIXlist is to take a vcf, convert it to a matrix, find informative sites
#   Goal of getLociPopAF is to calculate population allele frequency by loci
# TODO:
#     add better error handling
#     Change these functions to one call as an S4 object with slots for SNPmatrix, GT, PopAf, and d

#-----------------------------------------------------
# Read VCF and Go to SNP Matrix of Informative Sites
#------------------------------------------------------
vcf2infSNPmatrix <- function(vcffile, biallelic=T, diploid=T) {
  # Parses out the Format Field in a VCF file 
  #
  # Args:
  #   files: File path for vcf
  #   biallelic: Must be true. By default have set as false so user must think/confirm it is bialleic before running
  #   diploid: Must have heterozygous calls
  # Returns:
  #   A matrix of segregatig site/loci GT values as coded by adegenet (0=homoref, 1=het, 2=homoalt)
  
  ############################
  ##### Error handling #######
  ############################
  # Check to see if vcfR, TidyVerse, and DT are loaded
  # following this stackoverflow https://stackoverflow.com/questions/4090169/elegant-way-to-check-for-missing-packages-and-install-them
  list.of.packages <- c("vcfR", "tidyverse", "adegenet")
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) install.packages(new.packages) # here it will install if not already installed
  
  if (diploid != TRUE || biallelic != TRUE) {
    stop("This vcf must be biallelic and have heterozygote calls (i.e. called with a ploidy of 2")
  }
  
  # Check input
  if(is.null(vcffile)){
    stop("Must input a vcf file.")
  }
  
  ##################################################
  ##### Read in VCF and convert to SNPMatrix #######
  ##################################################
  
  datvcf <- read.vcfR(file=vcffile, verbose=T) # read vcf
  snpmatrix <- extract.gt(datvcf, element='GT', as.numeric=F)
  snpmatrix[snpmatrix == "0/0"] <- 0
  snpmatrix[snpmatrix == "0/1"] <- 1
  snpmatrix[snpmatrix == "1/1"] <- 2
  snpmatrix <- apply(snpmatrix, 2, function(x){as.numeric(x)})
  uninformativeloci <- apply(snpmatrix, 1, function(x){
    sum(x, na.rm=T) != 0 }) # find segregating sites
  snpmatrix_inform <- snpmatrix[uninformativeloci,] # subset to only sites with some variation
  
  
  #############################################################
  ##### Attach Positions to the VCF GT Informative Loci #######
  #############################################################
  POS <- data.frame(datvcf@fix, stringsAsFactors = F)
  POS <- POS[uninformativeloci,]
  POS <- as.numeric(POS$POS) # taking char to num -- fine in R, maybe not in Cpp
  snpmatrix_inform <- cbind(POS, snpmatrix_inform)
  return(snpmatrix_inform)
  
}




#-----------------------------------------------
# Calculate Pop AF per Loci 
#----------------------------------------------
getLociPopAF <- function(snpmatrix_rc_samploci){
  # Parses out the Format Field in a VCF file 
  #
  # Args:
  #   snpmatrix_rc_samploci: A SNP matrix of GT calls  
  #
  # Returns:
  #   A 2xN matrix of population allele frequencies, p & q
  
  ############################
  ##### Error handling #######
  ############################
  list.of.packages <- c("tidyverse")
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) install.packages(new.packages) # here it will install if not already installed
  
  # Check input
  if(is.null(snpmatrix_rc_samploci)){
    stop("Must input a snpmatrix that was generated using the vcf2infSNPmatrix function.")
  }
  
  
  ###################################
  #####      CODE for AF     #######
  ##################################
  popAF <- matrix(rep(NA, 2*ncol(snpmatrix_rc_samploci)), nrow=2) # init the matrix
  
  popAF[1,] <- apply(snpmatrix_rc_samploci[2:nrow(snpmatrix_rc_samploci),], 2, function(x){(2*(length(which(x==0))) + length(which(x==1)))/(2*length(x))}) # since we know homozygous ref is 0, so this counts as 2 As, and then we count hets and then divide by 2*possible alleles. Doing this roundabout way because we aren't in HWE
  popAF[2,] <- 1-popAF[1,] # now we have q
  
  return(popAF)
}


#-----------------------------------------------
# Calculate Genetic Distance
#----------------------------------------------
getLocigeneticdist <- function(snpmatrix_rc_samploci){
  # Parses out the Format Field in a VCF file 
  #
  # Args:
  #   snpmatrix_rc_samploci: A SNP matrix of GT calls  
  #
  # Returns:
  #   A vector of the genetic distance between consecutive loci 
  
  ############################
  ##### Error handling #######
  ############################
  list.of.packages <- c("tidyverse")
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) install.packages(new.packages) # here it will install if not already installed
  
  #####################################################
  #####      CODE for Genetic Distance by Loci #######
  ####################################################
  locigendist <- matrix(rep(NA, ncol(snpmatrix_rc_samploci)), nrow=1) # init the matrix
  locigendist[1] <- NA # redundant but to be explicit
  
  for(i in 2:nrow(snpmatrix_rc_samploci)){
    locigendist[i] <- snpmatrix_rc_samploci[i,1] - snpmatrix_rc_samploci[i-1,1]
  }
  
  return(locigendist)
}


