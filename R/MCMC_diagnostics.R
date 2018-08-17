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
