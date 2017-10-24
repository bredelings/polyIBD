#' @title Forward Algorithm for polyIBD
#' @description .....
#' @param file
#' @export
#' 

#------------------------------------------------
# Forward Algorithm
#------------------------------------------------

# note - m is now a vector of two values
# transProbs must be matrix with k+1 rows and k+1 columns (eventually an array)
# note, call this function with a particular subset of EmmissionLookupTable...
#      ForwardAlg(Genot, trans, emmission[[m1]][[m2]], m)

ForwardAlg <- function(GenotypeCompare, transProbs, EmissionLookUpTable, f) {
  
  # get simple parameters from inputs
  k <- nrow(transProbs)-1 # k is the maximum number of genotypes potentially shared between two samples 
  n <- length(GenotypeCompare) # length of the pariwise comparison for figuring out which level of the emission probability table
  
  # initialise matrix
  frwrd <- matrix(NA,k+1,n) # matrix for the different IBD states
  
  # loop through all loci
  for (i in 1:n) {
    if (i==1) { # first hmm state
      for (j in 1:(k+1)) {
        frwrd[j,1] <- dbinom(j-1, size=k, p=f)*EmissionLookUpTable[[j]][GenotypeCompare[i], i]
      }
    } else { # not first hmm state
      for (j in 1:(k+1)) {
        frwrd[j,i] <- sum(frwrd[,i-1]*transProbs[,j])*EmissionLookUpTable[[j]][GenotypeCompare[i], i]
      }
    }
  }
  return(frwrd)
}

#' @title Backward Algorithm for polyIBD
#' @description .....
#' @param file
#' @export
#' 

#------------------------------------------------
# Backward Algorithm
#------------------------------------------------

BackwardAlg <- function(GenotypeCompare, transProbs, EmissionLookUpTable) {
  
  # get simple parameters from inputs
  k <- nrow(transProbs)-1 # k is the maximum number of genotypes potentially shared between two samples 
  n <- length(GenotypeCompare) # length of the pariwise comparison for figuring out which level of the emission probability table
  
  # initialize matrix
  bkwrd <- matrix(1,k+1,n) # matrix for the different IBD states
  
  # loop through loci backwards from penultimate
  for (i in (n-1):1) {
    for (j in 1:(k+1)) {
      tmp <- 0
      for (j2 in 1:(k+1)) {
        tmp <- tmp + transProbs[j,j2]*bkwrd[j2,i+1]*EmissionLookUpTable[[j2]][GenotypeCompare[i+1], 1]
      }
      bkwrd[j,i] <- tmp
    }
  }
  return(bkwrd)
}

