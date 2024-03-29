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
#' @description Produces a simple MCMC trace plot of the parameter \code{m1},
#' which represents the COI of the first sample.
#'
#' @param x an object of class \code{polyIBD}, as produced by the function
#' \code{polyIBD::runMCMC}
#'
#' @export


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
  do.call(plot_trace, c(list(x=unclass(x$iterations$m1)), args))
}

#------------------------------------------------
#' @title Trace plot of m2
#'
#' @description Produces a simple MCMC trace plot of the parameter \code{m2},
#' which represents the COI of the second sample.
#'
#' @param x an object of class \code{polyIBD}, as produced by the function
#' \code{polyIBD::runMCMC}
#'
#' @export


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
  do.call(plot_trace, c(list(x=unclass(x$iterations$m2)), args))
}

#------------------------------------------------
#' @title Trace plot of f posterior probability
#'
#' @description Produces a simple MCMC trace plot of the parameter \code{f}, which represents the posterior probability of identity by descent between two samples.
#'
#' @param x an object of class \code{polyIBD}, as produced by the function \code{polyIBD::runMCMC}
#'
#' @export


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
  do.call(plot_trace, c(list(x=unclass(x$iterations$f)), args))
}

#------------------------------------------------
#' @title Trace plot of k
#'
#' @description Produces a simple MCMC trace plot of the parameter \code{k},
#' which represents the number of generations separating the two lineages
#' (holding recombination rate constant).
#'
#' @param x an object of class \code{polyIBD}, as produced by the function
#' \code{polyIBD::runMCMC}
#'
#' @export


plot_k <- function(x, ...) {

  # only works on objects of class polyIBD
  stopifnot(is.polyIBD(x))

  # get input arguments
  args <- list(...)
  argNames <- names(args)

  # set defaults on undefined arguments
  if (! "ylab" %in% argNames) {args$ylab <- "k"}
  if (! "main" %in% argNames) {args$main <- "k trace"}
  if (! "ylim" %in% argNames) {args$ylim <- c(0, max(unclass(x$iterations$k)))}

  # produce plot
  do.call(plot_trace, c(list(x = unclass(x$iterations$k)), args))
}

#------------------------------------------------
#' @title Plot marginal IBD matrix
#'
#' @description Plots the full posterior IBD matrix between two samples. SNPs
#' are arranged in order on the x-axis, and the IBD level is given on the y-axis. This represents the number of genotypes that are shared between samples, with the lowest level (IBD=0) representing zero shared genotypes and hence no IBD. The intensity of colour represents the probability of each IBD level at each locus.
#'
#' @param x an object of class \code{polyIBD}, as produced by the function
#' \code{polyIBD::runMCMC}
#' @param trueIBD option to overlay a line corresponding to the true IBD
#' (for example if using simulated data)
#'
#' @export


plot_IBD <- function(x, trueIBD = NULL, ...) {

  # only works on objects of class polyIBD
  stopifnot(is.polyIBD(x))

  # get input arguments
  args <- list(...)
  argNames <- names(args)

  # set defaults on undefined arguments
  if (! "xlab" %in% argNames) {args$xlab <- "Genomic position"}
  if (! "ylab" %in% argNames) {args$ylab <- "IBD level"}
  if (! "main" %in% argNames) {args$main <- "posterior mean IBD level"}
  if (! "bigplot" %in% argNames) {args$bigplot <- c(0.1, 0.85, 0.2, 0.8)}
  if (! "smallplot" %in% argNames) {args$smallplot <- c(0.9, 0.92, 0.3, 0.7)}
  if (! "legend.args" %in% argNames) {args$legend.args <- list(text  ="probability", side = 3, line  =0.5, adj  =0, cex  =0.8)}

  # get IBD matrix
  CHROM <- x$summary$IBD_marginal[,1]
  POS <- x$summary$IBD_marginal[,2]
  IBD <- as.matrix(x$summary$IBD_marginal[,-(1:2)])

  # produce plot
  do.call(IBDraster, c(list(x = IBD, CHROM = CHROM, POS = POS, whiteLine = trueIBD), args))
}

#------------------------------------------------
#' @title Default multi-panel plot for polyIBD class
#'
#' @description The default \code{plot()} function can be called on objects
#' returned from \code{polyIBD::runMCMC()} (which are of class \code{polyIBD}),
#' to produce a multi-panel plot giving important MCMC output.
#'
#' @param x an object of class \code{polyIBD}, as produced by the function
#' \code{polyIBD::runMCMC}
#'
#' @export

plot.polyIBD <- function(x, y = NULL, ...) {

  # temporarily switch to multi-panel
  layout( cbind(c(1,1,3,3,5,5,5), c(2,2,4,4,5,5,5)) )
  oldPar <- par(mar = c(4,4,2,1)) # default par
  on.exit({
    par(oldPar)
    layout(1)
  })

  # plot m1 and m2 trace
  plot_m1(x, ylim = c(0,max(x$iterations$m1) + 1))
  plot_m2(x, ylim = c(0,max(x$iterations$m2) + 1))

  # plot f
  plot_f(x)

  # plot k
  plot_k(x)

  # plot IBD matrix
  plot_IBD(x, trueIBD = y, ...)
}

#------------------------------------------------
# faster
# Produce multiple adjacent image plots separated by contig, with an associated
# legend. Makes use of functions from the fields pacakge, and hence includes many
# arguments also found in fields::image.plot.
# (not exported)

IBDraster <- function (x, col = NULL, breaks = NULL, nlevel = 64,
                       layout_mat = NULL, horizontal = FALSE,
                       legend.shrink = 0.9, legend.width = 1.2,
                       legend.mar = ifelse(horizontal, 3.1, 5.1),
                       legend.lab = NULL, legend.line = 2,  bigplot = NULL,
                       smallplot = NULL, lab.breaks = NULL, axis.args = NULL,
                       legend.args = NULL, legend.cex = 1, CHROM = NULL,
                       POS = NULL, whiteLine = NULL, ...) {

  #-----------------
  # set params

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
  layout_simple <- matrix(1:(mfrow[1]*mfrow[2]), mfrow[1], mfrow[2], byrow = TRUE)

  # set layout to simple grid if not specified
  if (is.null(layout_mat)) {
    layout_mat <- layout_simple
  }

  # set breaks evenly over data range if not specified
  info <- fields::imagePlotInfo(x, breaks = breaks, nlevel = nlevel)
  breaks <- info$breaks

  # get plotting limits
  temp <- fields::imageplot.setup(add = FALSE, legend.shrink = legend.shrink,
                          legend.width = legend.width, legend.mar = legend.mar,
                          horizontal = horizontal, bigplot = bigplot,
                          smallplot = smallplot)
  smallplot <- temp$smallplot
  bigplot <- temp$bigplot


  #-----------------
  # add legend

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


  #-----------------
  # add main image plots

  # extract basic genomic position parameters
  tab1 <- table(CHROM)
  nc <- length(tab1)
  cnames <- names(tab1)
  n <- as.vector(tab1)

  # main image plot
  par(old.par)
  old.par <- par(new = TRUE, plt = bigplot)
  image(1:sum(n), 0:(ncol(x) - 1), x, breaks = breaks, add = FALSE, col = col, axes=FALSE, ...)

  # add dotted white lines between contigs
  abline(v = cumsum(n)[-length(n)], col = "white", lty = 3)

  # add axes
  axis(1, at = cumsum(n) - n/2, labels = cnames)
  axis(2)

  # add white lines
  if (!is.null(whiteLine)) {
    x0 <- 0
    for (i in 1:nc) {
      lines((x0 + 1):(x0+n[i]), whiteLine[[i]], col="white", lwd=2)
      x0 <- x0 + n[i]
    }
  }


  #-----------------
  # tidy up

  # get position of current and next plot
  mfg <- par()$mfg
  thisPlot <- layout_mat[mfg[1], mfg[2]]
  if (max(layout_simple) > 1) {
    nextPlot <- which(layout_simple==(thisPlot+1), arr.ind=TRUE)[1,]
  }

  # switch back to old pars
  par(old.par)

  # set position of next plot
  if (max(layout_simple) > 1) {
    par(mfg=c(nextPlot[1], nextPlot[2], nrow(layout_mat), ncol(layout_mat)))
  }
  par(new=FALSE)
}



