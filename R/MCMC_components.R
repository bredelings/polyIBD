
#------------------------------------------------
# Trans Probs
#------------------------------------------------
# ------------------------------------------------------------------
#' @title Get transition probabilities
#'
#' @description Takes values of f, rho, and zmax. Produces rate matrix and calculates eigen values and vectors.
#'
#' @param f TODO
#' @export

getTransProbs <- function(f, rho, k, z_max) {

  # generate rate matrix
  z0 <- z_max
  z1 <- z_max + 1
  rateMat <- matrix(0, z1, z1)
  rateMat[cbind(1:z0, 1:z0 + 1)] <- (z0:1)*rho*k*f # fill in matrix with flows up from current state to next state
  rateMat[cbind(1:z0 + 1, 1:z0)] <- (1:z0)*rho*k*(1 - f) # ...
  rateMat[cbind(1:z1, 1:z1)] <- -rowSums(rateMat)

  # obtain Eigen values and vectors
  E <- eigen(t(rateMat))
  Esolve <- solve(E$vectors)

  return(
    list(
      Evalues =  E$values,
      Evectors = mat_to_Rcpp(E$vectors),
      Esolve =   mat_to_Rcpp(Esolve)
    )
  )
}



#------------------------------------------------
# Convergence Diagnostics
#------------------------------------------------

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
  geweke_z <- coda::geweke.diag(chain, frac1 = nburnin/(nburnin+nsamples), frac2 = nsamples/(nburnin+nsamples))$z

  if(is.na(geweke_z)){
    stop("NaN p-value was calculated from Geweke statistic")
  }

  # convert to p-value
  geweke_p <- 2*pnorm(abs(geweke_z), lower.tail=FALSE)

  # report convergence
  if (geweke_p > 0.05) {
    cat(paste0("convergence reached within defined burn-in period (Geweke p=", round(geweke_p, 3), ")"))
  } else {
    cat(paste0("WARNING: convergence not reached within defined burn-in period (Geweke p=", round(geweke_p,3), ")"))
  }

}



# TODO add more diagnostics
# https://www2.stat.duke.edu/courses/Fall09/sta290/Lectures/Diagnostics/param-diag.pdf
# consider adding ACF plots
# effectiveSize(theta.MCMC)
