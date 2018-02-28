
# ------------------------------------------------------------------
#' @title Get transition probabilities
#'
#' @description Takes values of f, rho, and zmax. Produces rate matrix and calculates eigen values and vectors.
#'
#' @param f TODO
#' @export

getTransProbs <- function(f, rho, zmax) {
    
    # define alpha from f and rho
    alpha <- rho*f/(1-f)
    
    # generate rate matrix
    z0 <- zmax
    z1 <- zmax+1
    rateMat <- matrix(0, z1, z1)
    rateMat[cbind(1:z0, 1:z0 + 1)] <- (z0:1)*alpha
    rateMat[cbind(1:z0 + 1, 1:z0)] <- (1:z0)*rho
    rateMat[cbind(1:z1, 1:z1)] <- -rowSums(rateMat)
    
    # obtain Eigen values and vectors
    E <- eigen(t(rateMat))
    Esolve <- solve(E$vectors)
    
    return(list(Evalues=E$values, Evectors=mat_to_Rcpp(E$vectors), Esolve=mat_to_Rcpp(Esolve)))
}

# -----------------------------------
# checkConvergence
# calculates Geweke statistic from a series of burn-in and sampling draws. Report whether burn-in length was sufficient based on this statistic.
# (not exported)

checkConvergence <- function(burnin, samples) {
  
  # get number of burnin and sampling iterations
  nburnin <- length(burnin)
  nsamples <- length(samples)
  
  # calculate Geweke diagnostic on combined chain
  chain <- coda::mcmc(c(burnin, samples))
  geweke_z <- coda::geweke.diag(chain, frac1=nburnin/(nburnin+nsamples), frac2=nsamples/(nburnin+nsamples))$z
  
  # convert to p-value
  geweke_p <- 2*pnorm(abs(geweke_z), lower.tail=FALSE)
  
  # report convergence
  if (geweke_p > 0.05) {
    cat(paste0("convergence reached within defined burn-in period (Geweke p=", round(geweke_p, 3), ")"))
  } else {
    cat(paste0("WARNING: convergence not reached within defined burn-in period (Geweke p=", round(geweke_p,3), ")"))
  }
  
}

#------------------------------------------------
# function for determining if object is of class polyIBD
#' @export

is.polyIBD <- function(x) {
  inherits(x, "polyIBD")
}

# -----------------------------------
# mat_to_Rcpp
# takes matrix as input, converts to list format for use within Rcpp code
# (not exported)

mat_to_Rcpp <- function(x) {
    return(split(x,f=1:nrow(x)))
}

# -----------------------------------
# Rcpp_to_mat
# Takes list format returned from Rcpp and converts to matrix.
# (not exported)

Rcpp_to_mat <- function(x) {
    ret <- matrix(unlist(x), nrow=length(x), byrow=TRUE)
    return(ret)
}
