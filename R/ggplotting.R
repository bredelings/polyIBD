#------------------------------------------------
#' @title polyibdggtheme
#'

polyibdggtheme <- ggplot2::theme(panel.background = ggplot2::element_rect(fill = "white", colour = "grey50"),
               panel.grid.major = ggplot2::element_blank(),
               panel.grid.minor = ggplot2::element_blank(),
               plot.title =       ggplot2::element_text(size = 14, face = "bold", hjust = 0.5, family = "Arial"),
               axis.title.x =     ggplot2::element_text(size=12, face="bold", family = "Arial"),
               axis.text.x =      ggplot2::element_text(size=10, family = "Arial", angle = 45, hjust = 1),
               axis.title.y =     ggplot2::element_text(size=12, face="bold", family = "Arial"),
               axis.text.y =      ggplot2::element_text(size=12, family = "Arial"),
               axis.ticks =       ggplot2::element_blank()
)

#------------------------------------------------
#' @title ggtrace
#' @importFrom ggplot2 ggplot aes


ggplot_trace <- function(data=NULL, mapping=NULL, stat="identity", position="identity", na.rm=FALSE, inherit.aes = TRUE, ...) {
    ggplot2::layer(data=data,
                   mapping=mapping,
                   stat=stat,
                   position = position,
                   geom = ggplot2::GeomPoint,
                   inherit.aes = inherit.aes,
                   params = list( na.rm = na.rm,
                                  colour = "#252525",
                                  fill = "#969696",
                                  alpha=0.8))
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
#' @importFrom magrittr %>%
#' @importFrom ggplot2 ggplot aes
#'
#' @export


ggplot_m1 <- function(x) {

  # only works on objects of class polyIBD
  stopifnot(is.polyIBD(x))

  # set up for plot
  x <- x %>%
    polyIBD::polyIBD_iterations2tidy(.) %>%
    dplyr::filter(param == "m1")

  xbreaks <- as.numeric(floor(quantile(x = x$iteration, probs = c(seq(0, 1, by = 0.2)))))
  ybreaks <- 1:(max(x$value) + 1)

  # call plot
  x %>%
    ggplot2::ggplot(mapping = aes(x=iteration, y=value)) + ggplot_trace() +
    ggplot2::scale_y_continuous(name = "M1 Estimate", breaks = ybreaks, limits = c(0, max(ybreaks))) +
    ggplot2::scale_x_continuous(name = "MCMC Chain Iteration", breaks = xbreaks) +
    ggplot2::ggtitle("M1 Trace") +
    polyibdggtheme
}


#------------------------------------------------
#' @title Trace ggplot of m2
#'
#' @description Produces a simple MCMC trace ggplot2::geom_plot of the parameter \code{m2}, which represents the COI of the second sample.
#'
#' @param x an object of class \code{polyIBD}, as produced by the function \code{polyIBD::runMCMC}
#'
#' @importFrom magrittr %>%
#' @importFrom ggplot2 ggplot aes
#'
#' @export

ggplot_m2 <- function(x) {

  # only works on objects of class polyIBD
  stopifnot(is.polyIBD(x))

  # set up for plot
  x <- x %>%
    polyIBD::polyIBD_iterations2tidy(.) %>%
    dplyr::filter(param == "m2")

  xbreaks <- as.numeric(floor(quantile(x = x$iteration, probs = c(seq(0, 1, by = 0.2)))))
  ybreaks <- 1:(max(x$value) + 1)

  # call plot
  x %>%
    ggplot2::ggplot(mapping = aes(x=iteration, y=value)) + ggplot_trace() +
    ggplot2::scale_y_continuous(name = "M2 Estimate", breaks = ybreaks, limits = c(0, max(ybreaks))) +
    ggplot2::scale_x_continuous(name = "MCMC Chain Iteration", breaks = xbreaks) +
    ggplot2::ggtitle("M2 Trace") +
    polyibdggtheme
}





#------------------------------------------------
#' @title Trace ggplot of f posterior prob
#'
#' @description Produces a simple MCMC trace ggplot2::geom_plot of the parameter \code{f}, which represents the posterior probability of identity by descent between two samples.
#'
#' @param x an object of class \code{polyIBD}, as produced by the function \code{polyIBD::runMCMC}
#'
#' @importFrom magrittr %>%
#' @importFrom ggplot2 ggplot aes
#'
#' @export


ggplot_f_postprob <- function(x) {

    # only works on objects of class polyIBD
    stopifnot(is.polyIBD(x))

    # set up for plot
    x <- x %>%
      polyIBD::polyIBD_iterations2tidy(.) %>%
      dplyr::filter(param == "f")

    xbreaks <- as.numeric(floor(quantile(x = x$iteration, probs = c(seq(0, 1, by = 0.2)))))

    # call plot
    x %>%
      ggplot2::ggplot(mapping = aes(x=iteration, y=value)) + ggplot_trace() +
      ggplot2::scale_y_continuous(name = "F (Pop) Estimate", limits = c(0,1)) +
      ggplot2::scale_x_continuous(name = "MCMC Chain Iteration", breaks = xbreaks) +
      ggplot2::ggtitle("F (Pop) Trace") +
      polyibdggtheme
  }


#------------------------------------------------
#' @title Trace ggplot of f individual
#'
#' @description Produces a simple MCMC trace ggplot2::geom_plot of the parameter \code{f}, which represents the individual draw of the posterior probability of identity by descent between two samples.
#'
#' @param x an object of class \code{polyIBD}, as produced by the function \code{polyIBD::runMCMC}
#'
#' @importFrom magrittr %>%
#'
#' @importFrom ggplot2 ggplot aes
#'
#' @export


ggplot_f <- function(x) {

  # only works on objects of class polyIBD
  stopifnot(is.polyIBD(x))

  # set up for plot
  x <- x %>%
    polyIBD::polyIBD_iterations2tidy(.) %>%
    dplyr::filter(param == "f_ind")

  xbreaks <- as.numeric(floor(quantile(x = x$iteration, probs = c(seq(0, 1, by = 0.2)))))


  # call plot
  x %>%
    ggplot2::ggplot(mapping = aes(x=iteration, y=value)) + ggplot_trace() +
    ggplot2::scale_y_continuous(name = "F (Ind) Estimate", limits = c(0,1)) +
    ggplot2::scale_x_continuous(name = "MCMC Chain Iteration", breaks = xbreaks) +
    ggplot2::ggtitle("F (Ind) Trace") +
    polyibdggtheme
}


#------------------------------------------------
#' @title Trace plot of k
#'
#' @description Produces a simple MCMC trace plot of the parameter \code{k}, which represents the number of generations separating the two lineages (holding the recombination rate constant).
#'
#' @param x an object of class \code{polyIBD}, as produced by the function \code{polyIBD::runMCMC}
#'
#' @importFrom magrittr %>%
#'
#' @importFrom ggplot2 ggplot aes
#'
#' @export


ggplot_k <- function(x) {

  # only works on objects of class polyIBD
  stopifnot(is.polyIBD(x))

  # set up for plot
  x <- x %>%
    polyIBD::polyIBD_iterations2tidy(.) %>%
    dplyr::filter(param == "k")

  xbreaks <- as.numeric(floor(quantile(x = x$iteration, probs = c(seq(0, 1, by = 0.2)))))
  ybreaks <- 1:(max(x$value) + 1)

  # call plot
  x %>%
    ggplot2::ggplot(mapping = aes(x=iteration, y=value)) + ggplot_trace() +
    ggplot2::scale_y_continuous(name = "K Estimate", breaks = ybreaks, limits = c(0, max(ybreaks))) +
    ggplot2::scale_x_continuous(name = "MCMC Chain Iteration", breaks = xbreaks) +
    ggplot2::ggtitle("K Trace") +
    polyibdggtheme
}



#------------------------------------------------
#' @title Plot marginal IBD matrix
#'
#' @description Plots the full posterior IBD matrix between two samples. SNPs are arranged in order on the x-axis, and the IBD level is given on the y-axis. This represents the number of genotypes that are shared between samples, with the lowest level (IBD=0) representing zero shared genotypes and hence no IBD. The intensity of colour represents the probability of each IBD level at each locus.
#'
#' @param x an object of class \code{polyIBD}, as produced by the function \code{polyIBD::runMCMC}
#' @param trueIBD option to overlay a line corresponding to the true IBD (for example if using simulated data)
#'
#' @importFrom ggplot2 ggplot, aes
#'
#' @export



# TODO is lag right here? I don't like this...recombination blocks versus SNP absolute position issue
# TODO this should have a better return -- am going to need that better return for selection anyways

ggplot_IBD <- function(x, trueIBD=NULL, ...) {

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

  plotobj <- ggplot2::ggplot() +
    ggplot2::geom_rect(data = IBDdflong, mapping = aes(
      xmin = start, xmax = end, ymin = Znum - 0.49, ymax = Znum + 0.49, fill = Prob)) +
    viridis::scale_fill_viridis("IBD Probability", option = "plasma", limits = c(0,1)) +
    ggplot2::scale_y_continuous("Number of IBD Genotypes", breaks = seq(1:max(IBDdflong$Znum+1))-1) +
    ggplot2::xlab("POS") +
    ggplot2::facet_grid(~CHROM) +
    polyibdggtheme

  if(!is.null(trueIBD)){
    plotobj <- plotobj + ggplot2::geom_line(data = trueIBD, aes(x = POS, y = z_true),
                                            colour = "#de2d26", size = 0.75)
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
#' @importFrom ggplot2 ggplot aes
#'
#' @export


ggplot_IBDraster <- function(x, trueIBD = NULL, truem1 = NULL,
                             truem2 = NULL, truef = NULL, truek = NULL, truefpop = NULL) {
  if(!is.null(truem1)){
    plotm1 <- ggplot_m1(x)
    plotm1 <- plotm1 + ggplot2::geom_hline(yintercept = truem1, colour = "#de2d26", size = 1.1)
  } else {
    plotm1 <- ggplot_m1(x)
  } # end m1 plot

  if(!is.null(truem2)){
    plotm2 <- ggplot_m2(x)
    plotm2 <- plotm2 + ggplot2::geom_hline(yintercept = truem2,
                                           colour = "#de2d26", size = 1.1)
  } else {
    plotm2 <- ggplot_m2(x)
  } # end m2 plot


  if(!is.null(truefpop)){
    plotfpop <- ggplot_f_postprob(x)
    plotfpop <- plotfpop + ggplot2::geom_hline(yintercept=truefpop,
                                               colour = "#de2d26", size = 1.1)
  } else {
    plotfpop <- ggplot_f_postprob(x)
  } # end fpop plot


  if(!is.null(truef)){
    plotf <- ggplot_f(x)
    plotf <- plotf + ggplot2::geom_hline(yintercept=truef,
                                         colour = "#de2d26", size = 1.1)
  } else {
    plotf <- ggplot_f(x)
  } # end f plot

  if(!is.null(truek)){
    plotkpop <- ggplot_k(x)
    plotkpop <- plotkpop + ggplot2::geom_hline(yintercept=truek,
                                               colour = "#de2d26", size = 1.1)
  } else {
    plotkpop <- ggplot_k(x)
  } # end k plot

  # placeholder for now
  temp <- tibble::tibble(x=1, y=1)
  plotk <- ggplot2::ggplot(data=temp, ggplot2::aes(x=x, y=y, label="placeholder")) +
    ggplot2::geom_text() + ggplot2::theme_minimal() + ggplot2::ggtitle("K (Ind) Trace")


  if(!is.null(trueIBD)){
    plotibd <- ggplot_IBD(x, trueIBD)
  } else {
    plotibd <- ggplot_IBD(x)
  } # end IBD plot

  lay <- rbind(c(1,1,2,2),
               c(3,3,4,4),
               c(5,5,6,6),
               c(7,7,7,7),
               c(7,7,7,7))

  gridExtra::grid.arrange(plotm1, plotm2, plotfpop, plotkpop, plotf, plotk, plotibd,
                          layout_matrix=lay)

}

