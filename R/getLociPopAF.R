#' @title polyIBD popAF
#' @description .....
#' @param file
#' @export

# TO DO
#   Need to make this object

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
  popAF <- matrix(rep(NA, 2*(ncol(snpmatrix_rc_samploci)-2)), nrow=2) # init the matrix, -2 for the chromosome and position
  
  popAF[1,] <- apply(snpmatrix_rc_samploci[,3:ncol(snpmatrix_rc_samploci)], 2, function(x){(2*(length(which(x==0))) + length(which(x==1)))/(2*length(x))}) # since we know homozygous ref is 0, so this counts as 2 As, and then we count hets and then divide by 2*possible alleles. Doing this roundabout way because we aren't in HWE
  popAF[2,] <- 1-popAF[1,] # now we have q
  
  return(popAF)
}
