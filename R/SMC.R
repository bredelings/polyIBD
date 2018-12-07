
#------------------------------------------------
#' @title Simulate ancestral recombination graph under continent-island model
#'
#' @description Simulate ancestral recombination graph under continent-island
#'   model. Assumes N haploid individuals each generation, with a probability m
#'   of originating from the continent each generation.
#'
#' @param n number of haploid samples.
#' @param loci genomic positions of loci.
#' @param N haploid population size.
#' @param rho recombination rate.
#' @param mu probability of migration from continent (aka mutation rate).
#' @param generations number of generations to simulate.
#'
#' @export

continentisland_ARG <- function(n = 2,
                             loci = seq(0,1e3,1e2),
                             N = 100,
                             rho = 1e-3,
                             mu = 1e-2,
                             generations = 1e2) {
  
  # checks on inputs
  assert_single_pos_int(n, zero_allowed = FALSE)
  assert_vector(loci)
  assert_int(loci)
  assert_single_pos_int(N, zero_allowed = FALSE)
  assert_single_pos(rho)
  assert_single_pos(mu, zero_allowed = FALSE)
  assert_bounded(mu, left = 0, right = 1, inclusive_left = FALSE, inclusive_right = FALSE)
  assert_single_pos_int(generations, zero_allowed = FALSE)
  
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

#' @param recom boolean matrix specifying whether there was a recombination
#'   between this locus and the previous locus.
#' @param recom_time matrix specifying recombination times.
#' @param event_time matrix of event timings.
#' @param coalescence boolean matrix specifying whether event is coalescence
#'   (TRUE) or migration (FALSE).
#' @param coalesce_with matrix specifying, in the case of coalescence, which
#'   other lineage to coalesce to.
#' @param migrate_with matrix specifying, in the case of migration, which member
#'   of the continent population is the original ancestor.
#' @param loci vector of genomic positions.
#' @param N haploid population size.
#' @param rho recombination rate.
#' @param m probability of migration from the continent per generation.
#'
#' @export

continentisland_SMC_conditional <- function(recom = NULL,
                                            recom_time = NULL,
                                            event_time = NULL,
                                            coalescence = NULL,
                                            coalesce_with = NULL,
                                            migrate_with = NULL,
                                            loci = seq(0,1e3,1e2),
                                            N = 1e2,
                                            rho = 1e-3,
                                            m = 1e-2) {
  
  # checks on inputs
  assert_vector(loci)
  assert_int(loci)
  assert_single_pos_int(N, zero_allowed = FALSE)
  assert_single_pos(rho)
  assert_single_pos(mu, zero_allowed = FALSE)
  assert_bounded(mu, left = 0, right = 1, inclusive_left = FALSE, inclusive_right = FALSE)
  
  # find whether this is the first draw from the conditional process
  first_draw <- is.null(event_time)
  
  # make copies of input matrices
  recom_input <- recom
  recom_time_input <- recom_time
  event_time_input <- event_time
  coalescence_input <- coalescence
  coalesce_with_input <- coalesce_with
  migrate_with_input <- migrate_with
  
  # process inputs if not first draw
  next_migrant <- 0
  if (!first_draw) {
    
    # get value of next migrant
    next_migrant <- max(migrate_with, na.rm = TRUE) + 1
    
    # replace NA with -1 before passing to C++
    recom_time_input[is.na(recom_time_input)] <- -1
    coalesce_with_input[is.na(coalesce_with_input)] <- -1
    migrate_with_input[is.na(migrate_with_input)] <- -1
    
    # get matrices into list format
    recom_input <- mat_to_Rcpp(recom_input)
    recom_time_input <- mat_to_Rcpp(recom_time_input)
    event_time_input <- mat_to_Rcpp(event_time_input)
    coalescence_input <- mat_to_Rcpp(coalescence_input)
    coalesce_with_input <- mat_to_Rcpp(coalesce_with_input)
    migrate_with_input <- mat_to_Rcpp(migrate_with_input)
  }
  
  # define argument list
  args <- list(first_draw = first_draw,
               recom = recom_input,
               recom_time = recom_time_input,
               event_time = event_time_input,
               coalescence = coalescence_input,
               coalesce_with = coalesce_with_input,
               migrate_with = migrate_with_input,
               loci = loci,
               N = N,
               rho = rho,
               m = m,
               next_migrant = next_migrant)
  
  # run efficient C++ code
  output_raw <- continentisland_SMC_conditional_cpp(args)
  
  # process output
  recom <- cbind(recom, output_raw$recom_new)
  recom_time <- cbind(recom_time, output_raw$recom_time_new)
  recom_time[recom == FALSE] <- NA
  event_time <- cbind(event_time, output_raw$event_time_new)
  coalescence <- cbind(coalescence, output_raw$coalescence_new)
  coalesce_with <- cbind(coalesce_with, output_raw$coalesce_with_new)
  coalesce_with[coalescence == FALSE] <- NA
  migrate_with <- cbind(migrate_with, output_raw$migrate_with_new)
  migrate_with[coalescence == TRUE] <- NA
  
  # compute ancestry matrix
  ancestry <- migrate_with
  n <- ncol(recom)
  for (j in 1:n) {
    for (i in 1:n) {
      y <- coalesce_with[,i]
      if (any(!is.na(y))) {
        ancestry[!is.na(y),i] <- ancestry[cbind(which(!is.na(y)), y[!is.na(y)]+1)]
      }
    }
  }
  
  # get ancestry spectrum matrix
  ancestry_spectrum <- t(apply(ancestry, 1, function(x) tabulate(match(x, unique(x)), nbins = n)))
  
  # get effective number of genotypes at every locus
  n_genotypes <- rowSums(ancestry_spectrum != 0)
  
  # return list
  ret <- list(recom = recom,
              recom_time = recom_time,
              event_time = event_time,
              coalescence = coalescence,
              coalesce_with = coalesce_with,
              migrate_with = migrate_with,
              ancestry = ancestry,
              ancestry_spectrum = ancestry_spectrum,
              n_genotypes = n_genotypes)
  return(ret)
}
