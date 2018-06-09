#------------------------------------------------
# The following commands ensure that package dependencies are listed in the NAMESPACE file.

#' @useDynLib polyIBD
#' @importFrom Rcpp evalCpp
#' @importFrom RcppParallel RcppParallelLibs
NULL

#------------------------------------------------
#' Run polyIBD MCMC using Rcpp functions
#'
#' @param dat TODO
#'
#' @export


runMCMC <- function(gtmat=NULL, p=NULL,  m_max=5,
                    k_max=50, rho=1e-5, 
                    burnin=1e2, samples=1e3, e1=0.05, e2=0.05, reportIteration=1e3) {
  
  # ------------------------------
  #         PROCESS INPUTS
  # ------------------------------
  require(vcfR)
  require(tidyverse)
  
  # -----------------------------------------------------
  # Read and check input
  #------------------------------------------------------
  #stopifnot(is.polyIBDinput(input))
  
  # TODO - input parameter checks
  # note - vcf must have 4 columns, samples in final two columns
  # better management of memory for p and vcf -- if you have it as one object, p continually gets copied and takes up way too much memory
  
  vcf <- gtmat
  p <- p
 
  
  # extract basic parameters
  tab1 <- table(vcf[,1]) # first column in this class is CHROM
  nc <- length(tab1)
  cnames <- names(tab1)
  n <- as.vector(tab1)
  
  # get distances between SNPs. Distance=-1 between contigs, indicating infinite distance
  SNP_dist <- diff(vcf[,2]) # second column in this class is POS
  SNP_dist[cumsum(n)[1:(nc-1)]] <- -1
  
  # compare two samples and save comparison type in vector x
  # x is an integer vector with values in 0:15. These values indicate genotype combinations that cycle through the four options: {missing, homo REF, het, homo ALT} in the first sample, then the same four options in the second sample, leading to 16 options in total
  x <- 4*(vcf[,3]+1) + (vcf[,4]+1)
  # fill in NA gaps
  x[is.na(vcf[,3]) & is.na(vcf[,4])] <- 0
  x[is.na(vcf[,3]) & vcf[,4] == 0] <- 1
  x[is.na(vcf[,3]) & vcf[,4] == 1] <- 2
  x[is.na(vcf[,3]) & vcf[,4] == 2] <- 3
  x[vcf[,3] == 0 & is.na(vcf[,4])] <- 4
  x[vcf[,3] == 1 & is.na(vcf[,4])] <- 8
  x[vcf[,3] == 2 & is.na(vcf[,4])] <- 12
  # define list of arguments to pass to Rcpp
  args <- list(x = x,
               p = unlist(p),
               rho = rho,
               SNP_dist = SNP_dist,
               m_max = m_max,
               k_max = k_max,
               burnin = burnin,
               samples = samples,
               e1 = e1,
               e2 = e2,
               reportIteration = reportIteration)
  
  # R functions to pass to Rcpp
  args_functions <- list(getTransProbs = polyIBD::getTransProbs)
  
  
  # ------------------------------
  #            RUN MCMC
  # ------------------------------
  
  # run efficient Rcpp function
  output_raw <- runMCMC_cpp(args, args_functions)
  
  # check for convergence
  checkConvergence(output_raw$logLike_burnin, output_raw$logLike)
  
  
  # ------------------------------
  #      SAVE MCMC RAW OUTPUT
  # ------------------------------
  
  # list of raw output
  raw_output <- list(logLike_burnin = coda::mcmc(output_raw$logLike_burnin),
                     logLike = coda::mcmc(output_raw$logLike),
                     m1 = coda::mcmc(output_raw$m1),
                     m2 = coda::mcmc(output_raw$m2),
                     f = coda::mcmc(output_raw$f),
                     f_ind = coda::mcmc(output_raw$f_ind),
                     k = coda::mcmc(output_raw$k),
                     #sim_trans_n = coda::mcmc(output_raw$sim_trans_n),
                     runTime = output_raw$runTime)
  
  
  # ------------------------------
  #   SAVE MCMC SUMMARY OUTPUT
  # ------------------------------
  
  # get marginal IBD matrix
  IBD_marginal <- t(Rcpp_to_mat(output_raw$IBD_marginal))
  colnames(IBD_marginal) <- paste0("z", 0:(ncol(IBD_marginal)-1))
  IBD_marginal <- cbind(vcf[,1:2], IBD_marginal)
  
  # get final acceptance rate
  accept_rate <- output_raw$accept_rate/samples
  
  # calculate quantiles over parameters
  quants <- t(mapply(function(x){quantile(x, probs=c(0.05, 0.5, 0.95))}, raw_output))
  quants <- quants[rownames(quants) %in% c("m1", "m2", "f", "f_ind", "k"),]
  
  # list of summary output
  summary_output <- list(IBD_marginal = IBD_marginal,
                         quantiles = quants,
                         accept_rate = accept_rate)
  
  
  # ------------------------------
  #            RETURN
  # ------------------------------
  
  # create output class polyIBD
  ret <- list(samples = colnames(vcf[3:4]),
              summary = summary_output,
              raw = raw_output)
  class(ret) <- "polyIBD"
  
  # return
  return(ret)
}

