#' @title polyIBD getLocigeneticdist
#' @description .....
#' @param file
#' @export

# TO DO
#   Need to make this object

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
  locigendist <- matrix(rep(NA, ncol(snpmatrix_rc_samploci)-2), nrow=1) # init the matrix, -2 for the chromosome and position columns
  locigendist[1] <- NA # redundant but to be explicit
  
  for(i in 2:nrow(snpmatrix_rc_samploci)){
    locigendist[i] <- snpmatrix_rc_samploci[i,2] - snpmatrix_rc_samploci[i-1,2] # column 2 contains chromosome position
  }
  
  return(locigendist)
}


