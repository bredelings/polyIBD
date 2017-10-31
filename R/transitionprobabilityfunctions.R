# ------------------------------------------------------------------
#' @title Calculate Transition Probabilities
#'
#' @description Take in model parameters and return transition probability matrix, evaluated at a given distance.
#'
#' @param rho rate of recombination
#' @param d genetic distance
#' @param f probability of IBD at equilibrium
#' @param zmax maximum level of IBD possible. For simple binary IBD/not-IBD model use zmax=1
#' @export

getTransProbs <- function(f, rho, zmax=1, d=1) {
  # define alpha from f and rho
  alpha <- rho*f/(1-f) 
  
  # generate rate matrix
  z0 <- zmax
  z1 <- zmax+1
  rateMat <- matrix(0, z1, z1)
  rateMat[cbind(1:z0, 1:z0 + 1)] <- (z0:1)*alpha
  rateMat[cbind(1:z0 + 1, 1:z0)] <- (1:z0)*rho
  rateMat[cbind(1:z1, 1:z1)] <- -rowSums(rateMat)
  
  # obtain Eigen values
  E <- eigen(t(rateMat))
  Esolve <- solve(E$vectors)
  
  # make transition matrix
  transProbs <- matrix(0, z1, z1)
  for (i in 1:z1) {
    # solve system
    constants <- E$vectors*outer(rep(1,z1),Esolve[,i])
    resMat <- constants*exp(outer(rep(1,z1),E$values)*d)
    transProbs[i,] <- rowSums(resMat)
  }
  
  return(transProbs)
}
