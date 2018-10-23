#------------------------------------------------
# The following commands ensure that package dependencies are listed in the NAMESPACE file.

#' @useDynLib polyIBD, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @importFrom RcppParallel RcppParallelLibs
NULL



#------------------------------------------------
#' Run polyIBD Stage I MCMC using Rcpp functions
#'
#' @param dat TODO
#'
#' @export


stgIrunMCMC <- function(sample = NULL,
                        input = NULL,
                        m_max = 5,
                        k_max = 10,
                        rho = 1e-5,
                        burnin = 1e2,
                        samples = 1e3,
                        e1 = 0.05, e2 = 0.05,
                        reportIteration = 1e3) {


  stopifnot(is.polyIBDinput(input))

  # TODO - input parameter checks
  # note - vcf must have 4 columns, samples in final two columns
  # better management of memory for p and vcf -- if you have it as one object, p continually gets copied and takes up way too much memory

  # ------------------------------
  #     Setup Input for Rcpp
  # ------------------------------

  Rcppcompat <- polyIBDinput_to_stgIrunMCMC_compat(polyIBDinput = input, sample=sample)

  # define list of arguments to pass to Rcpp
  args <- list(x = Rcppcompat[["x"]],
               p = unlist(Rcppcompat[["p"]]),
               rho = rho,
               SNP_dist = Rcppcompat[["SNP_dist"]],
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
  output_raw <- stgIrunMCMC_cpp(args, args_functions)

  # check for convergence
  checkConvergence(output_raw$logLike_burnin, output_raw$logLike)


  # ------------------------------
  #      SAVE MCMC RAW OUTPUT
  # ------------------------------

  # list of raw output
  raw_output <- list(logLike_burnin = coda::mcmc(output_raw$logLike_burnin),
                     logLike = coda::mcmc(output_raw$logLike),
                     m1 = coda::mcmc(output_raw$m1),
                     f = coda::mcmc(output_raw$f),
                     k = coda::mcmc(output_raw$k),
                     fws = coda::mcmc(output_raw$fws)
                     #sim_trans_n = coda::mcmc(output_raw$sim_trans_n),
  )


  # ------------------------------
  #   SAVE MCMC SUMMARY OUTPUT
  # ------------------------------

  # get marginal IBD matrix
  IBD_marginal <- t(Rcpp_to_mat(output_raw$IBD_marginal))
  colnames(IBD_marginal) <- paste0("z", 0:(ncol(IBD_marginal)-1))
  IBD_marginal <- cbind(input$CHROMPOS, IBD_marginal)

  # get effective MOI matrix
  effMOI <- t(Rcpp_to_mat(output_raw$effMOI))
  effMOI <- cbind(input$CHROMPOS, effMOI)


  # get final acceptance rate
  accept_rate <- output_raw$accept_rate/samples

  # get time it took MCMC to run
  runTime <- output_raw$runTime

  # calculate quantiles over parameters
  quants <- t(mapply(function(x){quantile(x, probs=c(0.025, 0.5, 0.975))}, raw_output))
  quants <- quants[rownames(quants) %in% c("m1", "f", "k", "fws"),]

  # list of summary output
  summary_output <- list(IBD_marginal = IBD_marginal,
                         effMOI = effMOI,
                         quantiles = quants,
                         accept_rate = accept_rate,
                         runTime = runTime)


  # ------------------------------
  #            RETURN
  # ------------------------------

  # create output class polyIBD
  ret <- list(samples = colnames(input$gtmatrix),
              summary = summary_output,
              iterations = raw_output)
  class(ret) <- "polyIBD"

  # return
  return(ret)
}




#------------------------------------------------
#' Run polyIBD Stage II MCMC using Rcpp functions
#'
#' @param dat TODO
#'
#' @export


runMCMC <- function(input = NULL,
                    retsmpl1, retsmpl2,
                    m_max = 5,
                    k_max = 10,
                    rho = 1e-5,
                    burnin = 1e2,
                    samples = 1e3,
                    e1 = 0.05, e2 = 0.05,
                    reportIteration = 1e3) {


  stopifnot(is.polyIBDinput(input))

  # TODO - input parameter checks
  # note - vcf must have 4 columns, samples in final two columns
  # better management of memory for p and vcf -- if you have it as one object, p continually gets copied and takes up way too much memory

  # ------------------------------
  #     Setup Input for Rcpp
  # ------------------------------

  Rcppcompat <- polyIBDinput_to_stgIIrunMCMC_compat(input, retsmpl1, retsmpl2)

  # define list of arguments to pass to Rcpp
  args <- list(x = Rcppcompat[["x"]],
               p = unlist(Rcppcompat[["p"]]),
               rho = rho,
               SNP_dist = Rcppcompat[["SNP_dist"]],
               m_max = m_max,
               k_max = k_max,
               burnin = burnin,
               samples = samples,
               e1 = e1,
               e2 = e2,
               reportIteration = reportIteration,
               m1_mat =  mat_to_Rcpp(t(Rcppcompat[["effm1_mat"]])), # transpose here so that loci are cols
               m2_mat =  mat_to_Rcpp(t(Rcppcompat[["effm2_mat"]]))
  )

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
                     k = coda::mcmc(output_raw$k)
                     #sim_trans_n = coda::mcmc(output_raw$sim_trans_n),
                     )


  # ------------------------------
  #   SAVE MCMC SUMMARY OUTPUT
  # ------------------------------

  # get marginal IBD matrix
  IBD_marginal <- t(Rcpp_to_mat(output_raw$IBD_marginal))
  colnames(IBD_marginal) <- paste0("z", 0:(ncol(IBD_marginal)-1))
  IBD_marginal <- cbind(input$CHROMPOS, IBD_marginal)

  # get final acceptance rate
  accept_rate <- output_raw$accept_rate/samples

  # get time it took MCMC to run
  runTime <- output_raw$runTime

  # calculate quantiles over parameters
  quants <- t(mapply(function(x){quantile(x, probs=c(0.025, 0.5, 0.975))}, raw_output))
  quants <- quants[rownames(quants) %in% c("m1", "m2", "f", "f_ind", "k"),]

  # list of summary output
  summary_output <- list(IBD_marginal = IBD_marginal,
                         quantiles = quants,
                         accept_rate = accept_rate,
                         runTime = runTime)


  # ------------------------------
  #            RETURN
  # ------------------------------

  # create output class polyIBD
  ret <- list(samples = colnames(input$gtmatrix),
              summary = summary_output,
              iterations = raw_output)
  class(ret) <- "polyIBD"

  # return
  return(ret)
}

