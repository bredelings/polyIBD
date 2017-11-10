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
#' runMCMC2()

runMCMC2 <- function(vcf, p, rho=1, m_max=5, burnin=1e2, samples=1e3) {
    
    # TODO - input parameter checks
    # note - vcf must have 4 columns, samples in final two columns
    
    # extract basic parameters
    L <- nrow(vcf)  # number of loci
    
    # compare two samples and save comparison type in vector x
    x <- 3*vcf[,3] + vcf[,4]
    
    # get distances between SNPs
    # TODO - infinite distances between chromosomes
    SNP_dist <- diff(vcf$POS)
    
    # define list of arguments
    args <- list(x=x,
                p=p,
                SNP_dist=SNP_dist,
                rho=rho,
                m_max=m_max,
                burnin=burnin,
                samples=samples)
    
    # run efficient Rcpp function
    output_raw <- runMCMC2_cpp(args)

    # post-processing of the raw output
    ret <- output_raw$x
    
    return(ret)
}

