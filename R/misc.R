
# ------------------------------------------------------------------
#' @title logit transform
#'
#' @description logit transform
#'
#' @param x value to be transformed
#' @export

logit <- function(x) {
    log(x) - log(1-x)
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
