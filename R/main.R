#------------------------------------------------
# The following commands ensure that package dependencies are listed in the NAMESPACE file.

#' @useDynLib polyIBD
#' @importFrom Rcpp evalCpp
#' @importFrom RcppParallel RcppParallelLibs
NULL

#------------------------------------------------
#' Run polyIBD MCMC using Rcpp functions
#'
#' Run polyIBD MCMC using Rcpp functions.
#'
#' @param dat TODO
#'
#' @export
#' @examples
#' runMCMC()

runMCMC <- function(vcf, p, m_max=5, rho_max=1e-5, burnin=1e2, samples=1e3, e1=0.05, e2=0.05, reportIteration=1e3) {
  
  # ------------------------------
  #         PROCESS INPUTS
  # ------------------------------
  
  # TODO - input parameter checks
  # note - vcf must have 4 columns, samples in final two columns
  
  # extract basic parameters
  L <- nrow(vcf)  # number of loci
  
  # compare two samples and save comparison type in vector x
  # x is an integer vector with values in 0:15. These values indicate genotype combinations that cycle through the four options: {missing, homo REF, het, homo ALT} in the first sample, then the same four options in the second sample, leading to 16 options in total
  x <- 4*(vcf[,3]+1) + (vcf[,4]+1)
  
  # get distances between SNPs
  # TODO - allow for multiple chromosomes, with infinite distances between chromosomes
  SNP_dist <- diff(vcf$POS)
  
  # define list of arguments to pass to Rcpp
  args <- list(x=x,
               p=p,
               rho_max=rho_max,
               SNP_dist=SNP_dist,
               m_max=m_max,
               burnin=burnin,
               samples=samples,
               e1=e1,
               e2=e2,
               reportIteration=reportIteration)
  
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
                     rho = coda::mcmc(output_raw$rho))
  
  
  # ------------------------------
  #   SAVE MCMC SUMMARY OUTPUT
  # ------------------------------
  
  # get marginal IBD matrix
  IBD_marginal <- Rcpp_to_mat(output_raw$IBD_marginal)
  colnames(IBD_marginal) <- rownames(vcf)
  rownames(IBD_marginal) <- paste0("z", 0:(nrow(IBD_marginal)-1))
  
  # list of summary output
  summary_output <- list(IBD_marginal = IBD_marginal)
  
  
  # ------------------------------
  #            RETURN
  # ------------------------------
  
  # create output class polyIBD
  ret <- list(summary = summary_output,
              raw = raw_output)
  class(ret) <- "polyIBD"
  
  # return
  return(ret)
}

