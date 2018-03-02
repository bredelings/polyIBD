
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
  
  # set defaults on undefined arguments
  if (! "xlab" %in% argNames) {args$xlab <- "SNP"}
  if (! "ylab" %in% argNames) {args$ylab <- "IBD level"}
  if (! "main" %in% argNames) {args$main <- "posterior IBD level"}
  
  # get IBD matrix
  IBD <- unclass(x$summary$IBD_marginal)
  
  # produce plot
  do.call(faster, c(list(x=1:ncol(IBD), y=0:(nrow(IBD)-1), z=t(IBD), bigplot=c(0.1, 0.85, 0.2, 0.8), smallplot=c(0.9, 0.92, 0.3, 0.7), legend.args=list(text="probability", side=3, line=0.5, adj=0, cex=0.8), whiteLine=trueIBD), args))
  
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
  plot_IBD(x, trueIBD=y, ...)
}

#------------------------------------------------
# faster
# Custom version of raster plot that is lightweight and doesn't suffer from problems relating to layout. Ordinary raster produces strange results when used as part of more complex layouts due to its method of setting and re-setting par(). This is avoided here by passing in an optional layout matrix that corresponds to the current layout. See ?fields::image.plot for details of what all arguments do.
# (not exported)

faster <- function (..., col=NULL, breaks=NULL, nlevel=64, layout_mat=NULL, horizontal=FALSE, legend.shrink=0.9, legend.width=1.2, legend.mar=ifelse(horizontal, 3.1, 5.1), legend.lab=NULL, legend.line=2,  bigplot=NULL, smallplot=NULL, lab.breaks=NULL, axis.args=NULL, legend.args=NULL, legend.cex=1, whiteLine=NULL) {
  
  # set defaults
  if (!is.null(breaks)) {
    nlevel <- length(breaks)-1
  }
  if (is.null(col)) {
    col <- viridis::plasma(nlevel)
  }
  nlevel <- length(col)
  
  if (is.null(legend.mar)) {
    legend.mar <- ifelse(horizontal, 3.1, 5.1)
  }
  old.par <- NULL
  
  # get layout matrix assuming simple increasing grid of values
  mfrow <- par()$mfrow
  layout_simple <- matrix(1:(mfrow[1]*mfrow[2]), mfrow[1], mfrow[2], byrow=TRUE)
  
  # set layout to simple grid if not specified
  if (is.null(layout_mat)) {
    layout_mat <- layout_simple
  }
  
  # set breaks evenly over data range if not specified
  nlevel <- length(col)
  info <- fields::imagePlotInfo(..., breaks = breaks, nlevel = nlevel)
  breaks <- info$breaks
  
  # get plotting limits
  temp <- fields::imageplot.setup(add = FALSE, legend.shrink = legend.shrink, 
                          legend.width = legend.width, legend.mar = legend.mar, 
                          horizontal = horizontal, bigplot = bigplot, smallplot = smallplot)
  smallplot <- temp$smallplot
  bigplot <- temp$bigplot
  
  # check legend will fit
  if ((smallplot[2] < smallplot[1]) | (smallplot[4] < smallplot[3])) {
    par(old.par)
    stop("plot region too small to add legend\n")
  }
  
  # make colour scale
  ix <- 1:2
  iy <- breaks
  nBreaks <- length(breaks)
  midpoints <- (breaks[1:(nBreaks - 1)] + breaks[2:nBreaks])/2
  iz <- matrix(midpoints, nrow = 1, ncol = length(midpoints))
  
  # add colour scale
  par(old.par)
  old.par <- par(pty = "m", plt = smallplot, err = -1)
  if (!horizontal) {
    image(ix, iy, iz, xaxt = "n", yaxt = "n", xlab = "", 
          ylab = "", col = col, breaks = breaks)
  }
  else {
    image(iy, ix, t(iz), xaxt = "n", yaxt = "n", xlab = "", 
          ylab = "", col = col, breaks = breaks)
  }
  
  # add numbers and box to scale
  if (!is.null(lab.breaks)) {
    axis.args <- c(list(side = ifelse(horizontal, 1, 4), 
                        mgp = c(3, 1, 0), las = ifelse(horizontal, 0, 2), 
                        at = breaks, labels = lab.breaks), axis.args)
  }
  else {
    axis.args <- c(list(side = ifelse(horizontal, 1, 4), 
                        mgp = c(3, 1, 0), las = ifelse(horizontal, 0, 2)), 
                   axis.args)
  }
  do.call("axis", axis.args)
  box()
  
  # add legend to scale
  if (!is.null(legend.lab)) {
    legend.args <- list(text = legend.lab, side = ifelse(horizontal, 
                                                         1, 4), line = legend.line, cex = legend.cex)
  }
  if (!is.null(legend.args)) {
    do.call(mtext, legend.args)
  }
  
  # main image plot
  old.par <- par(new = TRUE, plt = bigplot)
  image(..., breaks = breaks, add = FALSE, col = col)
  
  # add white line
  if (!is.null(whiteLine)) {
    lines(whiteLine, col="white", lwd=2)
  }
  
  # get position of current and next plot
  mfg <- par()$mfg
  thisPlot <- layout_mat[mfg[1], mfg[2]]
  if (max(layout_simple)>1) {
    nextPlot <- which(layout_simple==(thisPlot+1), arr.ind=TRUE)[1,]
  }
  
  # switch back to old pars
  par(old.par)
  
  # set position of next plot
  if (max(layout_simple)>1) {
    par(mfg=c(nextPlot[1], nextPlot[2], nrow(layout_mat), ncol(layout_mat)))
  }
  par(new=FALSE)
}

