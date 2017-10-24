#' @title polyIBD Transition Probabilities 
#' @description .....
#' @param rho
#' @param d
#' @param f
#' @param k
#' @export


#---------------------------------------
# Transition Probabilities
#---------------------------------------

getTransProbs <- function(f, rho, k, d) {
  # define alpha from f and rho
  alpha <- rho*f/(1-f) 
  
  # generate rate matrix
  rateMat <- matrix(0,k+1,k+1)
  rateMat[cbind(1:k,1:k+1)] <- (k:1)*alpha
  rateMat[cbind(1:k+1,1:k)] <- (1:k)*rho
  rateMat[cbind(1:(k+1),1:(k+1))] <- -rowSums(rateMat)
  
  # obtain Eigen values
  E <- eigen(t(rateMat))
  Esolve <- solve(E$vectors)
  
  # make transition matrix
  transProbs <- matrix(0,k+1,k+1)
  for (i in 1:(k+1)) {
    # solve system
    constants <- E$vectors*outer(rep(1,k+1),Esolve[,i])
    resMat <- constants*exp(outer(rep(1,k+1),E$values)*d)
    transProbs[i,] <- rowSums(resMat)
  }
  return(transProbs)
}
