
#------------------------------------------------
# plot_trace
# simple MCMC trace plot
# (not exported)

plot_trace <- function(x, ...) {
  
  # get input arguments
  args <- list(...)
  argNames <- names(args)
  
  # set defaults on undefined arguments
  if (! "pch" %in% argNames) {args$pch <- 20}
  if (! "col" %in% argNames) {args$col <- "#00000020"}
  if (! "xlab" %in% argNames) {args$xlab <- "iteration"}
  
  # produce plot
  do.call(plot, c(list(x=x), args))
}

#------------------------------------------------
#' @title Trace plot of m1
#'
#' @description Produces a simple MCMC trace plot of the parameter \code{m1}, which represents the COI of the first sample.
#'
#' @param x an object of class \code{polyIBD}, as produced by the function \code{polyIBD::runMCMC}
#'
#' @export
#' @examples

plot_m1 <- function(x, ...) {
  
  # only works on objects of class polyIBD
  stopifnot(is.polyIBD(x))
  
  # get input arguments
  args <- list(...)
  argNames <- names(args)
  
  # set defaults on undefined arguments
  if (! "ylab" %in% argNames) {args$ylab <- "m1"}
  if (! "main" %in% argNames) {args$main <- "m1"}
  
  # produce plot
  do.call(plot_trace, c(list(x=unclass(x$raw$m1)), args))
}

#------------------------------------------------
#' @title Trace plot of m2
#'
#' @description Produces a simple MCMC trace plot of the parameter \code{m2}, which represents the COI of the second sample.
#'
#' @param x an object of class \code{polyIBD}, as produced by the function \code{polyIBD::runMCMC}
#'
#' @export
#' @examples

plot_m2 <- function(x, ...) {
  
  # only works on objects of class polyIBD
  stopifnot(is.polyIBD(x))
  
  # get input arguments
  args <- list(...)
  argNames <- names(args)
  
  # set defaults on undefined arguments
  if (! "ylab" %in% argNames) {args$ylab <- "m2"}
  if (! "main" %in% argNames) {args$main <- "m2"}
  
  # produce plot
  do.call(plot_trace, c(list(x=unclass(x$raw$m2)), args))
}

#------------------------------------------------
#' @title Trace plot of f
#'
#' @description Produces a simple MCMC trace plot of the parameter \code{f}, which represents the mean probability of identity by descent between two samples.
#'
#' @param x an object of class \code{polyIBD}, as produced by the function \code{polyIBD::runMCMC}
#'
#' @export
#' @examples

plot_f <- function(x, ...) {
  
  # only works on objects of class polyIBD
  stopifnot(is.polyIBD(x))
  
  # get input arguments
  args <- list(...)
  argNames <- names(args)
  
  # set defaults on undefined arguments
  if (! "ylab" %in% argNames) {args$ylab <- "f"}
  if (! "main" %in% argNames) {args$main <- "f"}
  if (! "ylim" %in% argNames) {args$ylim <- c(0,1)}
  
  # produce plot
  do.call(plot_trace, c(list(x=unclass(x$raw$f)), args))
}

#------------------------------------------------
#' @title Trace plot of rho
#'
#' @description Produces a simple MCMC trace plot of the parameter \code{rho}, which represents the recombination rate.
#'
#' @param x an object of class \code{polyIBD}, as produced by the function \code{polyIBD::runMCMC}
#'
#' @export
#' @examples

plot_rho <- function(x, ...) {
  
  # only works on objects of class polyIBD
  stopifnot(is.polyIBD(x))
  
  # get input arguments
  args <- list(...)
  argNames <- names(args)
  
  # set defaults on undefined arguments
  if (! "ylab" %in% argNames) {args$ylab <- "rho"}
  if (! "main" %in% argNames) {args$main <- "rho"}
  if (! "ylim" %in% argNames) {args$ylim <- c(0, max(unclass(x$raw$rho)))}
  
  # produce plot
  do.call(plot_trace, c(list(x=unclass(x$raw$rho)), args))
}

#------------------------------------------------
#' Plot marginal IBD matrix
#'
#' text
#'
#' @param dat TODO
#'
#' @export
#' @examples

plot_IBD <- function(x, ...) {
  
  # only works on objects of class polyIBD
  stopifnot(is.polyIBD(x))
  
  # get input arguments
  args <- list(...)
  argNames <- names(args)
  
  # get IBD matrix
  IBD <- unclass(x$summary$IBD_marginal)
  
  # set defaults on undefined arguments
  if (! "x" %in% argNames) {args$x <- 1:ncol(IBD)}
  if (! "y" %in% argNames) {args$y <- 0:(nrow(IBD)-1)}
  if (! "col" %in% argNames) {args$col <- viridis::plasma(100)}
  if (! "xlab" %in% argNames) {args$xlab <- "SNP"}
  if (! "ylab" %in% argNames) {args$xlab <- "IBD level"}
  
  # produce plot
  do.call(image, c(list(z=t(IBD)), args))
  
  
  #image(1:ncol(IBD), 0:(nrow(IBD)-1), t(IBD), col=viridis::plasma(100), xlab="SNP", ylab="IBD level")
  
}
