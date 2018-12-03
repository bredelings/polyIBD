
#------------------------------------------------
#' @title Simulate ancestral recombination graph under Wright-Fisher model
#'
#' @description Simulate ancestral recombination graph under haploid
#'   Wright-Fisher model.
#'
#' @param n number of haploid samples.
#' @param loci genomic positions of loci.
#' @param N haploid population size.
#' @param rho recombination rate.
#' @param generations number of generations to simulate. If -1 then run until
#'   all loci have coalesced (or until hard limit of 1e9 generation is reached).
#'
#' @export

wrightfisher_ARG <- function(n = 2,
                             loci = seq(0,1000,100),
                             N = 100,
                             rho = 1e-3,
                             generations = -1) {
  
  # TODO - checks on inputs
  
  # define argument list
  args <- list(n = n,
               loci = loci,
               N = N,
               rho = rho,
               generations = generations)
  
  # run efficient C++ code
  output_raw <- wrightfisher_ARG_cpp(args)
  
  # process output
  coalesced <- Rcpp_to_mat(output_raw$coalesced)
  coalesced_with <- Rcpp_to_mat(output_raw$coalesced_with)
  coalesced_time <- Rcpp_to_mat(output_raw$coalesced_time) + 1
  coalesced_time[coalesced == 0] <- NA
  
  # return list
  ret <- list(coalesced = coalesced,
              coalesced_with = coalesced_with,
              coalesced_time = coalesced_time)
  return(ret)
}

#------------------------------------------------
#' @title Simulate from sequentially Markov coalescent under Wright-Fisher model
#'
#' @description Simulate from sequentially Markov coalescent approximation to
#'   the full ancestral recombination graph under haploid Wright-Fisher model.
#'   Note that because time is implemented as discrete this approximation will
#'   break down for large sample sizes.
#'
#' @param n number of haploid samples.
#' @param loci genomic positions of loci.
#' @param N haploid population size.
#' @param rho recombination rate.
#'
#' @export

wrightfisher_SMC <- function(n = 2,
                             loci = seq(0,1000,100),
                             N = 100,
                             rho = 1e-3) {
  
  # TODO - checks on inputs
  
  # define argument list
  args <- list(n = n,
               loci = loci,
               N = N,
               rho = rho)
  
  # run efficient C++ code
  output_raw <- wrightfisher_SMC_cpp(args)
  
  # process output
  coalesced_time <- t(Rcpp_to_mat(output_raw$coalesced_time)) + 1
  
  # return list
  ret <- list(coalesced_time = coalesced_time)
  return(ret)
}

#------------------------------------------------
#' @title Simulate from SMC under Wright-Fisher model by looping through
#'   discrete generations
#'
#' @description Simulate from sequentially Markov coalescent (SMC) approximation
#'   to the full ancestral recombination graph (ARG) under haploid Wright-Fisher
#'   model. Simulate using naive method of looping through discrete generations.
#'
#' @param n number of haploid samples.
#' @param loci genomic positions of loci.
#' @param N haploid population size.
#' @param rho recombination rate.
#' @param generations number of generations to simulate. If -1 then run until
#'   all loci have coalesced (or until hard limit of 1e9 generation is reached).
#'
#' @export

wrightfisher_SMC_naive <- function(n = 2,
                                   loci = seq(0,1000,100),
                                   N = 100,
                                   rho = 1e-3,
                                   generations = 10) {
  
  # TODO - checks on inputs
  
  # define argument list
  args <- list(n = n,
               loci = loci,
               N = N,
               rho = rho,
               generations = generations)
  
  # run efficient C++ code
  output_raw <- wrightfisher_SMC_naive_cpp(args)
  
  # process output
  coalesced <- Rcpp_to_mat(output_raw$coalesced)
  coalesced_time <- Rcpp_to_mat(output_raw$coalesced_time) + 1
  coalesced_time[coalesced == 0] <- NA
  
  # return list
  ret <- list(coalesced = coalesced,
              coalesced_time = coalesced_time)
  return(ret)
}

#------------------------------------------------
#' @title Simulate from SMC under Wright-Fisher model by layering haploid
#'   genotypes
#'
#' @description Simulate from sequentially Markov coalescent (SMC) approximation
#'   to the full ancestral recombination graph (ARG) under haploid Wright-Fisher
#'   model. Simulate by layering haploid genotypes one at a time.
#'
#' @param n number of haploid samples.
#' @param loci genomic positions of loci.
#' @param N haploid population size.
#' @param rho recombination rate.
#'
#' @export

wrightfisher_SMC_layer <- function(n = 2,
                                   loci = seq(0,1000,100),
                                   N = 100,
                                   rho = 1e-3,
                                   generations = NULL) {
  
  # TODO - checks on inputs
  
  # define argument list
  args <- list(n = n,
               loci = loci,
               N = N,
               rho = rho)
  
  # run efficient C++ code
  output_raw <- wrightfisher_SMC_layer_cpp(args)
  
  # process output
  coalesced_with <- t(Rcpp_to_mat(output_raw$coalesced_with))
  coalesced_with <- rbind(0, coalesced_with)
  coalesced_time <- t(Rcpp_to_mat(output_raw$coalesced_time) + 1)
  coalesced_time <- rbind(apply(coalesced_time, 2, max), coalesced_time)
  
  # return list
  ret <- list(coalesced_with = coalesced_with,
              coalesced_time = coalesced_time)
  return(ret)
}