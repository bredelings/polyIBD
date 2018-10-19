#------------------------------------------------
# S3 object functions
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
  cat(paste("burn-in iterations:\t", length(x$iterations$logLike_burnin)) ,"\n")
  cat(paste("sampling iterations:\t", length(x$iterations$logLike)) ,"\n")
  cat(paste("acceptance rate:\t", x$summary$accept_rate) ,"\n")
  cat(paste("run-time (seconds):\t", round(x$summary$runTime, 3)) ,"\n")
  cat("\n")

  # print posterior parameter summary
  cat("# Posterior estimates\n")
  quants <- round(x$summary$quantiles, 3)
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
  cat("-------------------------------------- \n")
  cat(paste("There are", ncol(x$gtmatrix), "Samples"), "\n")
  cat(paste(nrow(x$gtmatrix), "Biallelic SNPs"), "\n")
  cat("-------------------------------------- \n")
  # return invisibly
  invisible(x)
}


#------------------------------------------------
# overload summary() function to print summary only
#' @export
summary.polyIBDinput <- function(x, ...) {

  # print this output line
  cat("-------------------------------------- \n")
  cat(paste("There are", ncol(x$gtmatrix), "Samples"), "\n")
  cat(paste(nrow(x$gtmatrix), "Biallelic SNPs"), "\n")
  cat("-------------------------------------- \n")
  # return invisibly
  invisible(x)
}




#------------------------------------------------
# R <> Cpp compatibility
#------------------------------------------------

# -----------------------------------
# mat_to_Rcpp
# takes matrix as input, converts to list format for use within Rcpp code
# needed for trans prob eigen method
# (not exported)

mat_to_Rcpp <- function(x) {
  return(split(x,f=1:nrow(x)))
}

# -----------------------------------
# Rcpp_to_mat
# Takes list format returned from Rcpp and converts to matrix.
# needed for trans prob eigen method
# (not exported)

Rcpp_to_mat <- function(x) {
  ret <- matrix(unlist(x), nrow=length(x), byrow=TRUE)
  return(ret)
}

# -----------------------------------
# polyIBDinput_to_stgIrunMCMC_compat
# Takes the polyIBDinput and makes it easier to be parsed for the Rcpp args for the first stage MCMC.
# (not exported)
# TODO fix this input


polyIBDinput_to_stgIrunMCMC_compat <- function(polyIBDinput){

  p <- polyIBDinput[["p"]]
  gtmatrix <- polyIBDinput[["gtmatrix"]]

  # extract basic parameters
  CHROMtab <- table(polyIBDinput$CHROMPOS[, 1])
  nc <- length(CHROMtab)
  cnames <- names(CHROMtab)
  n <- as.vector(CHROMtab)

  # get distances between SNPs. Distance=-1 between contigs, indicating infinite distance
  SNP_dist <- diff(polyIBDinput$CHROMPOS[, 2]) # second column in this class is POS
  SNP_dist[cumsum(n)[1:(nc-1)]] <- -1


  # compare within sample and save comparison type in vector x
  # x is an integer vector with values in 0:3 These values indicate genotype combinations that cycle through the four options: {missing, homo REF, het, homo ALT}
  # 0 = {NA} - remember NA is read in as -1
  # 1 = {A}
  # 2 = {Aa}
  # 3 = {a}


  x <- gtmatrix[,1]+1

## return
ret <- list(p=p,
            x=x,
            SNP_dist = SNP_dist)

return(ret)

}


# -----------------------------------
# polyIBDinput_to_stgIIrunMCMC_compat
# Takes the polyIBDinput and makes it easier to be parsed for the Rcpp args for the second stage MCMC.
# (not exported)
# TODO fix this input

polyIBDinput_to_stgIIrunMCMC_compat <- function(polyIBDinput){
  p <- polyIBDinput[["p"]]
  gtmatrix <- polyIBDinput[["gtmatrix"]]

  # extract basic parameters
  CHROMtab <- table(polyIBDinput$CHROMPOS[, 1])
  nc <- length(CHROMtab)
  cnames <- names(CHROMtab)
  n <- as.vector(CHROMtab)

  # get distances between SNPs. Distance=-1 between contigs, indicating infinite distance
  SNP_dist <- diff(polyIBDinput$CHROMPOS[, 2]) # second column in this class is POS
  SNP_dist[cumsum(n)[1:(nc-1)]] <- -1

    # compare two samples and save comparison type in vector x
    # x is an integer vector with values in 0:15. These values indicate genotype combinations that cycle through the four options: {missing, homo REF, het, homo ALT} in the first sample, then the same four options in the second sample, leading to 16 options in total
    # 0 = {NA, NA}
    # 1 = {NA, A}
    # 2 = {NA, Aa}
    # 3 = {NA, a}
    # 4 = {A, NA}
    # 5 = {A,A}
    # 6 = {A,Aa}
    # 7 = {A,a}
    # 8 = {Aa, NA}
    # 9 = {Aa,A}
    # 10 = {Aa,Aa}
    # 11 = {Aa,a}
    # 12 = {a, NA}
    # 13 = {a,A}
    # 14 = {a,Aa}
    # 15 = {a,a}

    x <- 4*(gtmatrix[,1]+1) + (gtmatrix[,2]+1)


  ## return
  ret <- list(p=p,
              x=x,
              SNP_dist = SNP_dist)

  return(ret)

}


