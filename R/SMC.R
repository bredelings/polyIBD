
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
#' @title Simulate from conditional SMC under Wright-Fisher model
#'
#' @description Simulate from conditional SMC under Wright-Fisher model.
#'
#' @param n number of haploid samples.
#' @param loci genomic positions of loci.
#' @param N haploid population size.
#' @param rho recombination rate.
#'
#' @export

wrightfisher_SMC_conditional <- function(n = 2,
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
  output_raw <- wrightfisher_SMC_conditional_cpp(args)
  
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

#------------------------------------------------
#' @title Simulate ancestral recombination graph under continent-island model
#'
#' @description Simulate ancestral recombination graph under haploid
#'   continent-island model, which is equivalent to assuming infinite-alleles
#'   mutation.
#'
#' @param n number of haploid samples.
#' @param loci genomic positions of loci.
#' @param N haploid population size.
#' @param rho recombination rate.
#' @param mu probability of migration from continent (aka mutation rate).
#' @param generations number of generations to simulate. If -1 then run until
#'   all loci have coalesced (or until hard limit of 1e9 generation is reached).
#'
#' @export

continentisland_ARG <- function(n = 2,
                             loci = seq(0,1000,100),
                             N = 100,
                             rho = 1e-3,
                             mu = 1e-2,
                             generations = -1) {
  
  # TODO - checks on inputs
  
  # define argument list
  args <- list(n = n,
               loci = loci,
               N = N,
               rho = rho,
               mu = mu,
               generations = generations)
  
  # run efficient C++ code
  output_raw <- continentisland_ARG_cpp(args)
  
  # process output
  coalesced <- Rcpp_to_mat(output_raw$coalesced)
  migrated <- Rcpp_to_mat(output_raw$migrated)
  coalesced_with <- Rcpp_to_mat(output_raw$coalesced_with) + 1
  coalesced_with[!coalesced | migrated] <- NA
  coalesced_time <- Rcpp_to_mat(output_raw$coalesced_time) + 1
  coalesced_time[!coalesced | migrated] <- NA
  migrated_with <- Rcpp_to_mat(output_raw$migrated_with)
  migrated_with[!migrated | coalesced] <- NA
  migrated_time <- Rcpp_to_mat(output_raw$migrated_time) + 1
  migrated_time[!migrated | coalesced] <- NA
  
  ancestry <- migrated_with
  for (j in 1:n) {
    for (i in 1:n) {
      x <- migrated_with[i,]
      y <- coalesced_with[i,]
      ancestry[i,!is.na(y)] <- ancestry[cbind(y[!is.na(y)],which(!is.na(y)))]
    }
  }
  
  ancestry <- matrix(match(ancestry, unique(as.vector(ancestry))), nrow(ancestry))
  
  # return list
  ret <- list(coalesced_with = coalesced_with,
              coalesced_time = coalesced_time,
              migrated_with = migrated_with,
              migrated_time = migrated_time,
              ancestry = ancestry)
  return(ret)
}

#------------------------------------------------
#' @title Simulate from conditional SMC under continent-island model
#'
#' @description Simulate from conditional SMC under continent-island model.
#'
#' @param n number of haploid samples.
#' @param loci genomic positions of loci.
#' @param N haploid population size.
#' @param rho recombination rate.
#' @param mu probability of migration from continent (aka mutation rate).
#' @param generations (dummy argument)
#'
#' @export

continentisland_SMC_conditional <- function(n = 2,
                                loci = seq(0,1000,100),
                                N = 100,
                                rho = 1e-3,
                                mu = 1e-2,
                                generations = -1) {
  
  # TODO - checks on inputs
  
  # define argument list
  args <- list(n = n,
               loci = loci,
               N = N,
               rho = rho,
               mu = mu,
               generations = generations)
  
  # run efficient C++ code
  output_raw <- continentisland_SMC_conditional_cpp(args)
  
  # process output
  event_time <- t(Rcpp_to_mat(output_raw$event_time)) + 1
  coalesced <- t(Rcpp_to_mat(output_raw$coalesced))
  coalesced_with <- t(Rcpp_to_mat(output_raw$coalesced_with)) + 1
  coalesced_with[!coalesced] <- NA
  migrated_with <- t(Rcpp_to_mat(output_raw$migrated_with)) + 1
  migrated_with[coalesced] <- NA
  
  ancestry <- migrated_with
  for (j in 1:1) {
    for (i in 1:n) {
      y <- coalesced_with[i,]
      if (any(!is.na(y))) {
        ancestry[i,!is.na(y)] <- ancestry[cbind(y[!is.na(y)],which(!is.na(y)))]
      }
    }
  }
  
  ret <- list(event_time = event_time,
              coalesced = coalesced,
              coalesced_with = coalesced_with,
              migrated_with = migrated_with,
              ancestry = ancestry)
  return(ret)
}
