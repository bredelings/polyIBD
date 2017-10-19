#------------------------------------------------
# The following commands ensure that package dependencies are listed in the NAMESPACE file.

#' @useDynLib polyIBD
#' @importFrom Rcpp evalCpp
#' @importFrom RcppParallel RcppParallelLibs
NULL

#------------------------------------------------
#' Dummy function
#'
#' Simple dummy function used to demonstrate linking between R and C++. Returns vector of draws from standard normal distribution.
#'
#' @param n number of random draws
#'
#' @export
#' @examples
#' dummy1()

dummy1 <- function(n=10) {

    cat("R code working\n")

    # NOTES
    # one of the nice things about linking together R and C++ is that we can play to the strengths of both coding languages. For example, C++ is very fast and efficient, but it's a pain to import data and do simple checks on format etc. So the general plan here will be to do pre-processing in R, then send a clean list of arguments to the Rcpp function which does the heavy lifting, then finally tidy things up again in R.

    # example: check that n is positive
    stopifnot(n>0)

    # define mean and standard deviation of normal distribution
    mu <- 5
    sigma <- 1

    # create a final set of arguments as a list, and pass these to the Rcpp function
    args <- list(n=n,
                mu=mu,
                sigma=sigma)
    output_raw <- dummy1_cpp(args)

    # optionally do some post-processing of the raw output
    ret <- output_raw$x

    return(ret)
}

