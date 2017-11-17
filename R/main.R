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

runMCMC <- function(vcf, p, rho=1, m_max=5, burnin=1e2, samples=1e3, e1=0.05, e2=0.05, reportIteration=1e3) {
    
    # ------------------------------
    #         PROCESS INPUTS
    # ------------------------------
    
    # TODO - input parameter checks
    # note - vcf must have 4 columns, samples in final two columns
    
    # extract basic parameters
    L <- nrow(vcf)  # number of loci
    
    # compare two samples and save comparison type in vector x
    x <- 3*vcf[,3] + vcf[,4]
    
    # get distances between SNPs
    # TODO - allow for multiple chromosomes, with infinite distances between chromosomes?
    SNP_dist <- diff(vcf$POS)
    
    # define list of arguments
    args <- list(x=x,
                p=p,
                SNP_dist=SNP_dist,
                rho=rho,
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
    conv <- c(output_raw$logLike_burnin, output_raw$logLike)
    conv_z <- coda::geweke.diag(coda::mcmc(conv), frac1=burnin/(burnin+samples), frac2=samples/(burnin+samples))$z
    conv_p <- 2*pnorm(abs(conv_z), lower.tail=FALSE)
    
    # report convergence
    if (conv_p>0.05) {
        cat(paste0("convergence reached within defined burn-in period (Geweke p=", round(conv_p,3), ")"))
    } else {
        cat(paste0("WARNING: convergence not reached within defined burn-in period (Geweke p=", round(conv_p,3), ")"))
    }
    
    # ------------------------------
    #      SAVE MCMC RAW OUTPUT
    # ------------------------------
    
    # save raw chains
    logLike_burnin_chain <- coda::mcmc(output_raw$logLike_burnin)
    logLike_chain <- coda::mcmc(output_raw$logLike)
    m1_chain <- coda::mcmc(output_raw$m1)
    m2_chain <- coda::mcmc(output_raw$m2)
    f_chain <- coda::mcmc(output_raw$f)
    rho_chain <- coda::mcmc(output_raw$rho)
    
    # list of raw output
    raw_output <- list(logLike_burnin = logLike_burnin_chain,
                        logLike = logLike_chain,
                        m1 = m1_chain,
                        m2 = m2_chain,
                        f = f_chain,
                        rho = rho_chain)
    
    
    # ------------------------------
    #   SAVE MCMC SUMMARY OUTPUT
    # ------------------------------
    
    # get marginal IBD matrix
    IBD_marginal <- Rcpp_to_mat(output_raw$IBD_marginal)
    colnames(IBD_marginal) <- rownames(vcf)
    rownames(IBD_marginal) <- paste0("z",0:(nrow(IBD_marginal)-1))
    
    # list of summary output
    summary_output <- list(IBD_marginal = IBD_marginal)
    
    
    # ------------------------------
    #            RETURN
    # ------------------------------
    
    ret <- list(summary = summary_output,
                raw = raw_output)
    class(ret) <- "polyIBD"
    
    return(ret)
}

