#------------------------------------------------
# S4 object functions
#------------------------------------------------

#------------------------------------------------
# function for determining if object is of class polyIBD
#' @export

is.polyIBD <- function(x) {
  inherits(x, "polyIBD")
}

#------------------------------------------------
# overload print() function to print summary only
#' @export
print.polyIBD <- function(x, ...) {

  # print summary only
  summary(x)

  # return invisibly
  invisible(x)
}

#------------------------------------------------
# overload summary() function.
#' @export
summary.polyIBD <- function(x, ...) {

  # print MCMC summary
  cat("# MCMC summary\n")
  cat(paste("burn-in iterations:\t", length(x$raw$logLike_burnin)) ,"\n")
  cat(paste("sampling iterations:\t", length(x$raw$logLike)) ,"\n")
  cat(paste("acceptance rate:\t", x$summary$accept_rate) ,"\n")
  cat(paste("run-time (seconds):\t", round(x$raw$runTime, 3)) ,"\n")
  cat("\n")

  # print posterior parameter summary
  cat("# Posterior estimates\n")
  quants <- x$summary$quantiles
  print(quants)
}

#------------------------------------------------
# function for determining if object is of class polyIBD
#' @export

is.polyIBDinput <- function(x) {
  inherits(x, "polyIBDinput")
}


#------------------------------------------------
# overload print() function to print summary only
#' @export
print.polyIBDinput <- function(x, ...) {

  # print this output line
  cat("**************** \n")
  cat(paste("There are", length(x$samples), "Sample Combinations in this VCFsnpmatrix Object"), "\n")
  cat(paste(nrow(x$snpmatrix), "Variants per Sample Combination"), "\n")
  cat("**************** \n")
  # return invisibly
  invisible(x)
}


#------------------------------------------------
# overload summary() function to print summary only
#' @export
summary.polyIBDinput <- function(x, ...) {

  # print this output line
  cat("**************** \n")
  cat(paste("There are", length(x$samples), "Sample Combinations in this VCFsnpmatrix Object"), "\n")
  cat(paste(nrow(x$snpmatrix), "Variants per Sample Combination"), "\n")
  cat("**************** \n")
  # return invisibly
  invisible(x)
}


#------------------------------------------------
# R to Cpp compatibility
#------------------------------------------------

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
