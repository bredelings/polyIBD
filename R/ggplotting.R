#------------------------------------------------
# ggplot_trace
# simple MCMC trace plot
# (not exported)


ggplot_trace <- function(colour = "#252525", fill = "#969696", alpha = 0.8) {
    geom_point(colour = colour, fill="#969696", alpha=0.8) +
    theme(panel.background=element_rect(fill = "white", colour = "grey50"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.title = element_text(size = 14, face = "bold", hjust = 0.5, family = "Arial"),
          axis.title.x=element_text(size=12, face="bold", family = "Arial"),
          axis.text.x=element_text(size=10, family = "Arial", angle = 45, hjust = 1),
          axis.title.y=element_text(size=16.5, face="bold", family = "Arial"),
          axis.text.y=element_text(size=16, family = "Arial"),
          axis.ticks = element_blank()
    )


}



#------------------------------------------------
#' @title Trace ggplot of m1
#'
#' @description Produces a simple MCMC trace ggplot2::geom_plot of the parameter
#' \code{m1}, which represents the COI of the first sample.
#'
#' @param x an object of class \code{polyIBD}, as produced by the function
#' \code{polyIBD::runMCMC}
#'
#' @export


ggplot_m1 <- function(x, ...) {

  require(tidyverse)

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


ggplot_m2 <- function(x, ...) {

  require(tidyverse)

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
#' @title Trace ggplot of f posterior prob
#'
#' @description Produces a simple MCMC trace ggplot2::geom_plot of the parameter \code{f}, which represents the posterior probability of identity by descent between two samples.
#'
#' @param x an object of class \code{polyIBD}, as produced by the function \code{polyIBD::runMCMC}
#'
#' @export


ggplot_f_postprob <- function(x, ...) {

  require(tidyverse)

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
#' @title Trace ggplot of f individual
#'
#' @description Produces a simple MCMC trace ggplot2::geom_plot of the parameter \code{f}, which represents the individual draw of the posterior probability of identity by descent between two samples.
#'
#' @param x an object of class \code{polyIBD}, as produced by the function \code{polyIBD::runMCMC}
#'
#' @export


ggplot_f <- function(x, ...) {

  require(tidyverse)

  # only works on objects of class polyIBD
  stopifnot(is.polyIBD(x))

  # get input arguments
  args <- list(...)

  # produce df for ggplot
  dat <- data.frame(iteration=1:length(x$raw$f_ind), mcmcchain=unclass(x$raw$f_ind))

  # produce plot
  ggplot_trace(data=dat, mapping = aes(x=iteration, y=mcmcchain)) +
    scale_y_continuous(name="F Estimate", limits=c(0,1)) +
    scale_x_continuous(name="MCMC Chain Iteration", breaks=c(as.numeric(floor(quantile(x=c(1:length(x$raw$f_ind)), probs=seq(0,1,0.25)))))) +
    ggtitle(label="F (Ind) Trace")
}



#------------------------------------------------
#' @title Trace plot of k
#'
#' @description Produces a simple MCMC trace plot of the parameter \code{k}, which represents the number of generations separating the two lineages (holding the recombination rate constant).
#'
#' @param x an object of class \code{polyIBD}, as produced by the function \code{polyIBD::runMCMC}
#'
#' @export


ggplot_k <- function(x, ...) {

  require(tidyverse)

  # only works on objects of class polyIBD
  stopifnot(is.polyIBD(x))

  # get input arguments
  args <- list(...)

  # produce df for ggplot
  dat <- data.frame(iteration=1:length(x$raw$k), mcmcchain=unclass(x$raw$k))

  # produce plot
  ggplot_trace(data=dat, mapping = aes(x=iteration, y=mcmcchain)) +
    scale_y_continuous(name="k Estimate") +
    scale_x_continuous(name="MCMC Chain Iteration", breaks=c(as.numeric(floor(quantile(x=c(1:length(x$raw$k)), probs=seq(0,1,0.25)))))) +
    ggtitle(label="k Trace")
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



# TODO is lag right here? I don't like this...recombination blocks versus SNP absolute position issue


ggplot_IBD <- function(x, trueIBD=NULL, ...) {

  require(tidyverse)
  require(viridis)

  # only works on objects of class polyIBD
  stopifnot(is.polyIBD(x))


  # get IBD matrix
  CHROM <- x$summary$IBD_marginal[,1]
  POS <- x$summary$IBD_marginal[,2]
  IBD <- as.matrix(x$summary$IBD_marginal[,-(1:2)])

  IBDdf <- cbind.data.frame(CHROM, POS, IBD)
  IBDdflong <- tidyr::gather(data=IBDdf, key="Z", value="Prob", 3:ncol(IBDdf))
  IBDdflong$Znum <- as.numeric(gsub("z", "", IBDdflong$Z))

  # split by chrom to avoid lag pos from diff chrom
  IBDdflonglist <- split(IBDdflong, f=IBDdflong$CHROM)
  IBDdflonglist <- lapply(IBDdflonglist, function(df){
    df$start <- dplyr::lag(df$POS)
    df$end <- df$POS
    return(df)
  })
  IBDdflong <- do.call("rbind", IBDdflonglist)

  # filter unneccessary Znumbers
  filtdf <- aggregate(IBDdflong$Prob, list(factor(IBDdflong$Znum)), sum)
  if(any(filtdf == 0)){
    filtdf <- filtdf[which(filtdf[,2] == 0), 1]
    IBDdflong <- IBDdflong %>% dplyr::filter(! Znum %in% filtdf )
  }

  plotobj <- ggplot() +
    geom_rect(data=IBDdflong, mapping=aes(xmin=start, xmax=end, ymin=Znum-0.49, ymax=Znum+0.49, fill=Prob)) +
    viridis::scale_fill_viridis("IBD Probability", alpha=0.2, option="plasma", limits=c(0,1)) +
    scale_y_continuous("Number of IBD Genotypes", breaks = seq(1:max(IBDdflong$Znum+1))-1) +
    xlab("POS") +
    facet_grid(~CHROM) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.ticks = element_blank(),
          axis.title.y = element_text(size=18, face="bold", family = "Arial"),
          axis.text.y = element_text(size=16, face="bold", family = "Arial", margin = margin(r = 1)),
          axis.title.x = element_text(size=12, face="bold", family = "Arial"),
          axis.text.x = element_text(size=9, family = "Arial", angle = 45, margin = margin(r = 1)),
          strip.text.x = element_text(size =12, face="bold", family = "Arial"),
          legend.text=element_text(size = 10.5, face="bold", family = "Arial"),
          legend.title=element_text(size = 11, face="bold", family = "Arial"))

  if(!is.null(trueIBD)){
    plotobj <- plotobj + geom_line(data=trueIBD, aes(x=POS, y=z_true), colour="#de2d26", size=0.75)
  }
  return(plotobj)
}


#------------------------------------------------
#' @title Plot polyIBD outputs collectively
#'
#' @description Plots the MCMC chains for M1, M2, f, and k (generations), as well as the posterior distribution for the IBD matrix.
#'
#' @param x an object of class \code{polyIBD}, as produced by the function \code{polyIBD::runMCMC}
#' @param trueIBD option to overlay a line corresponding to the true IBD (for example if using simulated data)
#'
#' @export

# TODO -- switch this to grid extra

ggplot_IBDraster <- function(x, trueIBD=NULL, truem1=NULL,
                             truem2=NULL, truef=NULL, truek=NULL) {
  if(!is.null(truem1)){
    plotm1 <- ggplot_m1(x)
    plotm1 <- plotm1 + geom_hline(yintercept=truem1, colour="#de2d26", size=1.1)
  } else {
    plotm1 <- ggplot_m1(x)
  } # end m1 plot

  if(!is.null(truem2)){
    plotm2 <- ggplot_m2(x)
    plotm2 <- plotm2 + geom_hline(yintercept=truem2, colour="#de2d26", size=1.1)
  } else {
    plotm2 <- ggplot_m2(x)
  } # end m2 plot

  if(!is.null(truef)){
    plotf <- ggplot_f(x)
    plotf <- plotf + geom_hline(yintercept=truef, colour="#de2d26", size=1.1)
  } else {
    plotf <- ggplot_f(x)
  } # end f plot

  if(!is.null(truek)){
    plotk <- ggplot_k(x)
    plotk <- plotk + geom_hline(yintercept=truek, colour="#de2d26", size=1.1)
  } else {
    plotk <- ggplot_k(x)
  } # end k plot

  if(!is.null(trueIBD)){
    plotibd <- ggplot_IBD(x, trueIBD)
  } else {
    plotibd <- ggplot_IBD(x)
  } # end IBD plot

  ggpubr::ggarrange(ggpubr::ggarrange(plotm1, plotm2,  plotf, plotk, nrow=2, ncol = 2),
                    plotibd,
                    nrow = 2
  )
}


#------------------------------------------------
#' @title Plot marginal IBD matrix with true IBD option and VCF that would have been generated for the two samples for comparison
#'
#' @description Plots the full posterior IBD matrix between two samples. SNPs are arranged in order on the x-axis, and the IBD level is given on the y-axis. This represents the number of genotypes that are shared between samples, with the lowest level (IBD=0) representing zero shared genotypes and hence no IBD. The intensity of colour represents the probability of each IBD level at each locus.
#'
#' @param x an object of class \code{polyIBD}, as produced by the function \code{polyIBD::runMCMC}
#' @param trueIBD option to overlay a line corresponding to the true IBD (for example if using simulated data)
#'
#' @export

ggplot_SNPs_IBD <- function(x, trueIBD, snps){

  library(tidyverse)
  library(grid)

  snps <- snps[[1]] # only need first element in vcf
  # get SNP matrix
  CHROM <- snps[,1]
  POS <- snps[,2]
  vars <- as.matrix(snps[,-(1:2)])
  varsdf <- cbind.data.frame(CHROM, POS, vars)

  varsdf$status <- ifelse(varsdf$Sample1 == varsdf$Sample2, "Concord", "Discord")
  varsdf$status[varsdf$Sample1 == 1 | varsdf$Sample2 == 1] <- "Het"

  # split by chrom to avoid lag pos from diff chrom
  varsdflist <- split(varsdf, f=varsdf$CHROM)
  varsdflist <- lapply(varsdflist, lagPOS <- function(df){
    df <- df %>%
      dplyr::mutate(start = lag(POS)) %>%
      dplyr::mutate(end = POS) %>%
      dplyr::mutate(statusfct = factor(status, levels=c("Concord", "Discord", "Het")))
    df$start[1] = df$end[1] # fix lag NA and just make it the first instiation

    return(df)
  }
  )

  varsdf <- do.call("rbind", varsdflist)

  ### now call IBD plotter
  plotIBD <- ggplot_IBD(x, trueIBD)

  plotsnps <- ggplot() +
    geom_rect(data=varsdf, mapping=aes(xmin=start, xmax=end, ymin=0.5, ymax=1.5, fill=statusfct)) +
    scale_fill_manual("Support", values=c("#005AC8", "#AA0A3C", "#F0F032", "#cccccc")) +
    facet_grid(~CHROM) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.ticks = element_blank(),
          axis.text.x = element_text(size=9, family = "Arial", angle = 45),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          strip.text.x = element_text(size =12, face="bold", family = "Arial"),
          legend.text=element_text(size = 10.5, face="bold", family = "Arial"),
          legend.title=element_text(size = 11, face="bold", family = "Arial")
    )

  grid::grid.newpage()
  grid::grid.draw(rbind(ggplotGrob(plotsnps), ggplotGrob(plotIBD), size = "last"))

}

