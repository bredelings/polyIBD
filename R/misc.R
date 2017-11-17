
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
