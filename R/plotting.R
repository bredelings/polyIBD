
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
  if (! "main" %in% argNames) {args$main <- "m1 trace"}
  
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
  if (! "main" %in% argNames) {args$main <- "m2 trace"}
  
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
  if (! "main" %in% argNames) {args$main <- "f trace"}
  if (! "ylim" %in% argNames) {args$ylim <- c(0,1)}
  
  # produce plot
  do.call(plot_trace, c(list(x=unclass(x$raw$f)), args))
}

#------------------------------------------------
#' @title Trace plot of rho
#'
#' @description Produces a simple MCMC trace plot of the parameter \code{rho}, which represents the inverse of the average length of a recombinant block, and is a function of both the recombination rate and the number of generations separating the two lineages.
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
  if (! "main" %in% argNames) {args$main <- "rho trace"}
  if (! "ylim" %in% argNames) {args$ylim <- c(0, max(unclass(x$raw$rho)))}
  
  # produce plot
  do.call(plot_trace, c(list(x=unclass(x$raw$rho)), args))
}

#------------------------------------------------
#' @title Plot marginal IBD matrix
#'
#' @description Plots the full posterior IBD matrix between two samples. SNPs are arranged in order on the x-axis, and the IBD level is given on the y-axis. This represents the number of genotypes that are shared between samples, with the lowest level (IBD=0) representing zero shared genotypes and hence no IBD. The intensity of colour represents the probability of each IBD level at each locus.
#'
#' @param x an object of class \code{polyIBD}, as produced by the function \code{polyIBD::runMCMC}
#' @param trueIBD option to overlay a line corresponding to the true IBD (for example if using simulated data)
#'
#' @export
#' @examples

plot_IBD <- function(x, trueIBD=NULL, ...) {
  
  # only works on objects of class polyIBD
  stopifnot(is.polyIBD(x))
  
  # get input arguments
  args <- list(...)
  argNames <- names(args)
  
  # get IBD matrix
  IBD <- unclass(x$summary$IBD_marginal)
  IBD <- IBD[nrow(IBD):1,] # flip vertically
  
  # produce raster image from IBD matrix
  r <- raster::raster(IBD, xmn=0.5, xmx=ncol(IBD)+0.5, ymn=-0.5, ymx=nrow(IBD)-0.5)
  
  # set defaults on undefined arguments
  if (! "col" %in% argNames) {args$col <- viridis::plasma(100)}
  if (! "xlab" %in% argNames) {args$xlab <- "SNP"}
  if (! "ylab" %in% argNames) {args$ylab <- "IBD level"}
  if (! "asp" %in% argNames) {args$asp <- NA}
  
  # produce plot
  do.call(raster::image, c(list(x=r), args))
  #do.call(raster::plot, c(list(x=r), args))
  
  # overlay true IBD
  if (!is.null(trueIBD)) {
    lines(trueIBD, col="white", lwd=2)
  }
}

#------------------------------------------------
#' @title Default multi-panel plot for polyIBD class
#'
#' @description The default \code{plot()} function can be called on objects returned from \code{polyIBD::runMCMC()} (which are of class \code{polyIBD}), to produce a multi-panel plot giving important MCMC output.
#'
#' @param x an object of class \code{polyIBD}, as produced by the function \code{polyIBD::runMCMC}
#'
#' @export

plot.polyIBD <- function(x, y=NULL, ...) {
  
  # temporarily switch to multi-panel
  layout( cbind(c(1,1,3,3,5,5,5), c(2,2,4,4,5,5,5)) )
  oldPar <- par(mar=c(4,4,2,1))
  on.exit({
    par(oldPar)
    layout(1)
  })
  
  # plot m1 and m2 trace
  plot_m1(x, ylim=c(0,max(x$raw$m1)+1))
  plot_m2(x, ylim=c(0,max(x$raw$m2)+1))
  
  # plot f
  plot_f(x)
  
  # plot rho
  plot_rho(x)
  
  # plot IBD matrix
  plot_IBD(x, trueIBD=y, main="posterior IBD level")
}

