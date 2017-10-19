#' @title polyIBD vcf2SNPmatrix
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


