##################################
###   Base Plots for PolyIBD   ###
##################################

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
  if (! "xlab" %in% argNames) {args$xlab <- "Genomic position"}
  if (! "ylab" %in% argNames) {args$ylab <- "IBD level"}
  if (! "main" %in% argNames) {args$main <- "posterior mean IBD level"}
  if (! "bigplot" %in% argNames) {args$bigplot <- c(0.1, 0.85, 0.2, 0.8)}
  if (! "smallplot" %in% argNames) {args$smallplot <- c(0.9, 0.92, 0.3, 0.7)}
  if (! "legend.args" %in% argNames) {args$legend.args <- list(text="probability", side=3, line=0.5, adj=0, cex=0.8)}
  
  # get IBD matrix
  CHROM <- x$summary$IBD_marginal[,1]
  POS <- x$summary$IBD_marginal[,2]
  IBD <- as.matrix(x$summary$IBD_marginal[,-(1:2)])
  
  # produce plot
  do.call(IBDraster, c(list(x=IBD, CHROM=CHROM, POS=POS, whiteLine=trueIBD), args))
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
# Produce multiple adjacent image plots separated by contig, with an associated legend. Makes use of functions from the fields pacakge, and hence includes many arguments also found in fields::image.plot.
# (not exported)

IBDraster <- function (x, col=NULL, breaks=NULL, nlevel=64, layout_mat=NULL, horizontal=FALSE, legend.shrink=0.9, legend.width=1.2, legend.mar=ifelse(horizontal, 3.1, 5.1), legend.lab=NULL, legend.line=2,  bigplot=NULL, smallplot=NULL, lab.breaks=NULL, axis.args=NULL, legend.args=NULL, legend.cex=1, CHROM=NULL, POS=NULL, whiteLine=NULL, ...) {
  
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
  layout_simple <- matrix(1:(mfrow[1]*mfrow[2]), mfrow[1], mfrow[2], byrow=TRUE)
  
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
                          horizontal = horizontal, bigplot = bigplot, smallplot = smallplot)
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
  image(1:sum(n), 0:(ncol(x)-1), x, breaks = breaks, add = FALSE, col = col, axes=FALSE, ...)
  
  # add dotted white lines between contigs
  abline(v=cumsum(n)[-length(n)], col="white", lty=3)
  
  # add axes
  axis(1, at=cumsum(n)-n/2, labels=cnames)
  axis(2)
  
  # add white lines
  if (!is.null(whiteLine)) {
    x0 <- 0
    for (i in 1:nc) {
      lines((x0+1):(x0+n[i]), whiteLine[[i]], col="white", lwd=2)
      x0 <- x0 + n[i]
    }
  }
  
  
  #-----------------
  # tidy up
  
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


##################################
###   GGplots for PolyIBD      ###
##################################
#------------------------------------------------
# ggplot_trace
# simple MCMC trace plot
# (not exported)


ggplot_trace <- function(mapping = NULL, data = NULL, stat = "identity", position = "identity",  inherit.aes = F, ...) {
  # dependencies
  list.of.packages <- c("tidyverse")
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) install.packages(new.packages) # here it will install if not already installed
  library("tidyverse")
  
  # produce plot
  ggplot(mapping = mapping,
         data = data,
         stat = stat,
         position = position, 
         inherit.aes = inherit.aes) + 
    geom_point(colour="#252525", fill="#969696", alpha=0.8) +
    theme(panel.background=element_rect(fill = "white", colour = "grey50"), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          plot.title = element_text(size = 14, face = "bold", hjust = 0.5, family = "Arial"),
          axis.title.x=element_text(size=12, face="bold", family = "Arial"), 
          axis.text.x=element_text(size=10, family = "Arial", angle = 45, hjust = 1), 
          axis.title.y=element_text(size=12, face="bold", family = "Arial"), 
          axis.text.y=element_text(size=10, family = "Arial"),
          axis.ticks = element_blank()
    )
  
  
}



#------------------------------------------------
#' @title Trace ggplot of m1
#'
#' @description Produces a simple MCMC trace ggplot2::geom_plot of the parameter \code{m1}, which represents the COI of the first sample.
#'
#' @param x an object of class \code{polyIBD}, as produced by the function \code{polyIBD::runMCMC}
#'
#' @export
#' @examples

ggplot_m1 <- function(x, ...) {
  
  # only works on objects of class polyIBD
  stopifnot(is.polyIBD(x))
  
  # get input arguments
  args <- list(...)
  
  # produce df for ggplot
  dat <- data.frame(iteration=1:length(x$raw$m1), mcmcchain=unclass(x$raw$m1))
  
  # produce plot
  ggplot_trace(data=dat, mapping = aes(x=iteration, y=mcmcchain)) + 
    scale_y_continuous(name="M1 Estimate", breaks=seq(1:max(unclass(x$raw$m1)+0.5)), 
                       limits=c(0.5,max(unclass(x$raw$m1)+0.5))) +
    scale_x_continuous(name="MCMC Chain Iteration", breaks=c(as.numeric(floor(quantile(x=c(1:length(x$raw$m1)), probs=seq(0,1,0.25)))))) +
    ggtitle(label="M1 Trace")
}

#------------------------------------------------
#' @title Trace ggplot of m2
#'
#' @description Produces a simple MCMC trace ggplot2::geom_plot of the parameter \code{m2}, which represents the COI of the second sample.
#'
#' @param x an object of class \code{polyIBD}, as produced by the function \code{polyIBD::runMCMC}
#'
#' @export
#' @examples

ggplot_m2 <- function(x, ...) {
  
  # only works on objects of class polyIBD
  stopifnot(is.polyIBD(x))
  
  # get input arguments
  args <- list(...)
  
  # produce df for ggplot
  dat <- data.frame(iteration=1:length(x$raw$m2), mcmcchain=unclass(x$raw$m2))
  
  # produce plot
  ggplot_trace(data=dat, mapping = aes(x=iteration, y=mcmcchain)) + 
    scale_y_continuous(name="M2 Estimate", breaks=seq(1:max(unclass(x$raw$m2)+0.5)), 
                       limits=c(0.5,max(unclass(x$raw$m2)+0.5))) +
    scale_x_continuous(name="MCMC Chain Iteration", breaks=c(as.numeric(floor(quantile(x=c(1:length(x$raw$m2)), probs=seq(0,1,0.25)))))) +
    ggtitle(label="M2 Trace")
}


#------------------------------------------------
#' @title Trace ggplot of f
#'
#' @description Produces a simple MCMC trace ggplot2::geom_plot of the parameter \code{f}, which represents the mean probability of identity by descent between two samples.
#'
#' @param x an object of class \code{polyIBD}, as produced by the function \code{polyIBD::runMCMC}
#'
#' @export
#' @examples

ggplot_f <- function(x, ...) {
  
  # only works on objects of class polyIBD
  stopifnot(is.polyIBD(x))
  
  # get input arguments
  args <- list(...)
  
  # produce df for ggplot
  dat <- data.frame(iteration=1:length(x$raw$f), mcmcchain=unclass(x$raw$f))
  
  # produce plot
  ggplot_trace(data=dat, mapping = aes(x=iteration, y=mcmcchain)) + 
    scale_y_continuous(name="F Estimate", limits=c(0,1)) +
    scale_x_continuous(name="MCMC Chain Iteration", breaks=c(as.numeric(floor(quantile(x=c(1:length(x$raw$f)), probs=seq(0,1,0.25)))))) +
    ggtitle(label="F Trace")
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

ggplot_rho <- function(x, ...) {
  
  # only works on objects of class polyIBD
  stopifnot(is.polyIBD(x))
  
  # get input arguments
  args <- list(...)
  
  # produce df for ggplot
  dat <- data.frame(iteration=1:length(x$raw$rho), mcmcchain=unclass(x$raw$rho))
  
  # produce plot
  ggplot_trace(data=dat, mapping = aes(x=iteration, y=mcmcchain)) + 
    scale_y_continuous(name="Rho Estimate") +
    scale_x_continuous(name="MCMC Chain Iteration", breaks=c(as.numeric(floor(quantile(x=c(1:length(x$raw$rho)), probs=seq(0,1,0.25)))))) +
    ggtitle(label="Rho Trace")
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


ggplot_IBD <- function(x, trueIBD=NULL, ...) {
  # only works on objects of class polyIBD
  stopifnot(is.polyIBD(x))
  
  list.of.packages <- c("tidyverse", "viridis")
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) install.packages(new.packages) # here it will install if not already installed
  lapply(list.of.packages, function(x){require(x, character.only = TRUE)})
  
  # get input arguments
  args <- list(...)
  
  # get IBD matrix
  CHROM <- x$summary$IBD_marginal[,1]
  POS <- x$summary$IBD_marginal[,2]
  IBD <- as.matrix(x$summary$IBD_marginal[,-(1:2)])
  
  IBDdf <- data.frame(cbind(CHROM, POS, IBD))
  IBDdflong <- tidyr::gather(data=IBDdf, key="Z", value="Prob", 3:ncol(IBDdf))
  IBDdflong$Znum <- as.numeric(gsub("z", "", IBDdflong$Z))
  
  ggplot() +
    geom_tile(data=IBDdflong, aes(x=factor(POS), y=Znum, fill=Prob)) +
    facet_grid(~ CHROM, scales="free_x", space="free") +
    scale_fill_continuous("", low="#a1d99b", high="#003c30") +
    geom_hline(yintercept=seq(0.5, (max(IBDdflong$Znum)-0.5), by=1), color="black", size=0.5) +
    scale_fill_viridis(alpha=0.2, option="plasma") +
    scale_y_continuous("Number of IBD Genotypes", breaks = seq(1:max(IBDdflong$Znum+1))-1) +
    xlab("POS") +
    guides(title="IBD Probability", labels = paste("0%", "100%")) + 
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          panel.background = element_blank(),
          axis.ticks = element_blank(),
          axis.text.x = element_blank(),
          #axis.text.x = element_text(size=9, family = "Arial", angle = 45),
          axis.title.y = element_text(size=14, face="bold", family = "Arial"),
          axis.title.x = element_text(size=12, face="bold", family = "Arial"))
  
  # add in an if statement to handle the trueIBD info
  # fix the x-axis 
  
}



