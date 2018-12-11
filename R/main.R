#------------------------------------------------
# The following commands ensure that package dependencies are listed in the NAMESPACE file.

#' @useDynLib polyIBD, .registration = TRUE
#' @importFrom Rcpp evalCpp
NULL

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
#' @param m probability of migration from continent (aka mutation rate).
#' @param generations number of generations to simulate.
#'
#' @export

continentisland_ARG <- function(n = 2,
                                loci = seq(0,1e3,1e2),
                                N = 100,
                                rho = 1e-3,
                                m = 1e-2,
                                generations = 1e2) {
  
  # checks on inputs
  assert_single_pos_int(n, zero_allowed = FALSE)
  assert_vector(loci)
  assert_int(loci)
  assert_single_pos_int(N, zero_allowed = FALSE)
  assert_single_pos(rho)
  assert_single_pos(m, zero_allowed = FALSE)
  assert_bounded(m, left = 0, right = 1, inclusive_left = FALSE, inclusive_right = FALSE)
  assert_single_pos_int(generations, zero_allowed = FALSE)
  
  # define argument list
  args <- list(n = n,
               loci = loci,
               N = N,
               rho = rho,
               m = m,
               generations = generations)
  
  # run efficient C++ code
  output_raw <- continentisland_ARG_cpp(args)
  
  # process output
  recom <- t(Rcpp_to_mat(output_raw$recom))
  recom_time <- t(Rcpp_to_mat(output_raw$recom_time))
  event_time <- t(Rcpp_to_mat(output_raw$event_time))
  active <- t(Rcpp_to_mat(output_raw$active))
  coalescence <- t(Rcpp_to_mat(output_raw$coalescence))
  coalesce_target <- t(Rcpp_to_mat(output_raw$coalesce_target))
  coalesce_target[active | !coalescence] <- NA
  migrate_target <- t(Rcpp_to_mat(output_raw$migrate_target))
  migrate_target[active | coalescence] <- NA
  
  # compute ancestry matrix
  ancestry <- migrate_target
  for (j in 1:n) {
    for (i in 1:n) {
      y <- coalesce_target[,i]+1
      if (any(!is.na(y))) {
        ancestry[!is.na(y),i] <- ancestry[cbind(which(!is.na(y)), y[!is.na(y)])]
      }
    }
  }
  
  # get numbers to be increasing over the matrix
  ancestry <- matrix(match(ancestry, unique(as.vector(ancestry))), nrow(ancestry))
  
  # get ancestry spectrum matrix
  ancestry_spectrum <- t(apply(ancestry, 1, function(x) tabulate(match(x, unique(x)), nbins = n)))
  
  # get effective number of genotypes at every locus
  n_genotypes <- rowSums(ancestry_spectrum != 0)
  
  # return list
  ret <- list(recom = recom,
              recom_time = recom_time,
              event_time = event_time,
              active = active,
              coalescence = coalescence,
              coalesce_target = coalesce_target,
              migrate_target = migrate_target,
              ancestry = ancestry,
              ancestry_spectrum = ancestry_spectrum,
              n_genotypes = n_genotypes)
  return(ret)
}


#------------------------------------------------
#' @title Simulate from conditional SMC under continent-island model
#'
#' @description Simulate from a conditional version of the SMC algorithm under
#'   continent-island model. This allows genotypes to be built up sequentially.
#'   The first six inputs should be left as \code{NULL} when producing the first
#'   genotype, and subsequently the output of the previous function call (which
#'   includes six objects with the same names) should be passed in.

#' @param recom boolean matrix specifying whether there was a recombination
#'   between this locus and the previous locus.
#' @param recom_time matrix specifying recombination times.
#' @param event_time matrix of event timings.
#' @param coalescence boolean matrix specifying whether event is coalescence
#'   (TRUE) or migration (FALSE).
#' @param coalesce_target matrix specifying, in the case of coalescence, which
#'   other lineage to coalesce to.
#' @param migrate_target matrix specifying, in the case of migration, which
#'   member of the continent population is the original ancestor.
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
                                            coalesce_target = NULL,
                                            migrate_target = NULL,
                                            loci = seq(0,1e3,1e2),
                                            N = 1e2,
                                            rho = 1e-3,
                                            m = 1e-2) {
  
  # checks on inputs
  assert_vector(loci)
  assert_int(loci)
  assert_single_pos_int(N, zero_allowed = FALSE)
  assert_single_pos(rho)
  assert_single_pos(m, zero_allowed = FALSE)
  assert_bounded(m, left = 0, right = 1, inclusive_left = FALSE, inclusive_right = FALSE)
  
  # find whether this is the first draw from the conditional process
  first_draw <- is.null(event_time)
  
  # make copies of input matrices
  recom_input <- recom
  recom_time_input <- recom_time
  event_time_input <- event_time
  coalescence_input <- coalescence
  coalesce_target_input <- coalesce_target
  migrate_target_input <- migrate_target
  
  # process inputs if not first draw
  next_migrant <- 0
  if (!first_draw) {
    
    # get value of next migrant
    next_migrant <- max(migrate_target, na.rm = TRUE) + 1
    
    # replace NA with -1 before passing to C++
    recom_time_input[is.na(recom_time_input)] <- -1
    coalesce_target_input[is.na(coalesce_target_input)] <- -1
    migrate_target_input[is.na(migrate_target_input)] <- -1
    
    # get matrices into list format
    recom_input <- mat_to_Rcpp(recom_input)
    recom_time_input <- mat_to_Rcpp(recom_time_input)
    event_time_input <- mat_to_Rcpp(event_time_input)
    coalescence_input <- mat_to_Rcpp(coalescence_input)
    coalesce_target_input <- mat_to_Rcpp(coalesce_target_input)
    migrate_target_input <- mat_to_Rcpp(migrate_target_input)
  }
  
  # define argument list
  args <- list(first_draw = first_draw,
               recom = recom_input,
               recom_time = recom_time_input,
               event_time = event_time_input,
               coalescence = coalescence_input,
               coalesce_target = coalesce_target_input,
               migrate_target = migrate_target_input,
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
  coalesce_target <- cbind(coalesce_target, output_raw$coalesce_target_new)
  coalesce_target[coalescence == FALSE] <- NA
  migrate_target <- cbind(migrate_target, output_raw$migrate_target_new)
  migrate_target[coalescence == TRUE] <- NA
  
  # compute ancestry matrix
  ancestry <- migrate_target
  n <- ncol(recom)
  for (j in 1:n) {
    for (i in 1:n) {
      y <- coalesce_target[,i]+1
      if (any(!is.na(y))) {
        ancestry[!is.na(y),i] <- ancestry[cbind(which(!is.na(y)), y[!is.na(y)])]
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
              coalesce_target = coalesce_target,
              migrate_target = migrate_target,
              ancestry = ancestry,
              ancestry_spectrum = ancestry_spectrum,
              n_genotypes = n_genotypes)
  return(ret)
}

#------------------------------------------------
#' @title Drop genotype from conditional SMC under continent-island model
#'
#' @description Drop genotype from conditional SMC under continent-island model.
#' 
#' @param i which genotype to drop.
#' @param recom boolean matrix specifying whether there was a recombination
#'   between this locus and the previous locus.
#' @param recom_time matrix specifying recombination times.
#' @param event_time matrix of event timings.
#' @param coalescence boolean matrix specifying whether event is coalescence
#'   (TRUE) or migration (FALSE).
#' @param coalesce_target matrix specifying, in the case of coalescence, which
#'   other lineage to coalesce to.
#' @param migrate_target matrix specifying, in the case of migration, which
#'   member of the continent population is the original ancestor.
#'
#' @export

continentisland_SMC_conditional_drop <- function(i = 1,
                                                 recom,
                                                 recom_time,
                                                 event_time,
                                                 coalescence,
                                                 coalesce_target,
                                                 migrate_target) {
  
  
  # TODO
  return(-9)
}

#------------------------------------------------
#' @title Inference from conditional SMC under continent-island model
#'
#' @description Inference from conditional SMC under continent-island model
#' 
#' @param x data in the form of a sequence of homo/het calls.
#' @param n complexity of infection of sample.
#' @param loci vector of genomic positions.
#' @param N haploid population size.
#' @param rho recombination rate.
#' @param m probability of migration from the continent per generation.
#'
#' @export

SMC_MCMC <- function(x,
                     n = 3,
                     loci = seq(0,1e3,1e2),
                     N = 1e2,
                     rho = 1e-3,
                     m = 1e-2) {
  
  
  # checks on inputs
  assert_vector(x)
  assert_in(x, c(0,0.5,1))
  assert_single_pos_int(n, zero_allowed = FALSE)
  assert_vector(loci)
  assert_int(loci)
  assert_same_length(x, loci)
  assert_single_pos_int(N, zero_allowed = FALSE)
  assert_single_pos(rho)
  assert_single_pos(m, zero_allowed = FALSE)
  assert_bounded(m, left = 0, right = 1, inclusive_left = FALSE, inclusive_right = FALSE)
  
  # generate initial sample configuration
  recom <- NULL
  recom_time <- NULL
  event_time <- NULL
  coalescence <- NULL
  coalesce_target <- NULL
  migrate_target <- NULL
  for (i in 1:n) {
    sample_init <- continentisland_SMC_conditional(recom = recom,
                                                   recom_time = recom_time,
                                                   event_time = event_time,
                                                   coalescence = coalescence,
                                                   coalesce_target = coalesce_target,
                                                   migrate_target = migrate_target,
                                                   loci = loci,
                                                   N = N,
                                                   rho = rho,
                                                   m = m)
    recom <- sample_init$recom
    recom_time <- sample_init$recom_time
    event_time <- sample_init$event_time
    coalescence <- sample_init$coalescence
    coalesce_target <- sample_init$coalesce_target
    migrate_target <- sample_init$migrate_target
  }
  
  # get value of next migrant
  next_migrant <- max(migrate_target, na.rm = TRUE) + 1
  
  # replace NA with -1 before passing to C++
  recom_time[is.na(recom_time)] <- -1
  coalesce_target[is.na(coalesce_target)] <- -1
  migrate_target[is.na(migrate_target)] <- -1
  
  # get matrices into list format
  recom <- mat_to_Rcpp(recom)
  recom_time <- mat_to_Rcpp(recom_time)
  event_time <- mat_to_Rcpp(event_time)
  coalescence <- mat_to_Rcpp(coalescence)
  coalesce_target <- mat_to_Rcpp(coalesce_target)
  migrate_target <- mat_to_Rcpp(migrate_target)
  
  # define argument list
  args <- list(x = match(x, c(0,0.5,1)),
               n = n,
               recom = recom,
               recom_time = recom_time,
               event_time = event_time,
               coalescence = coalescence,
               coalesce_target = coalesce_target,
               migrate_target = migrate_target,
               loci = loci,
               N = N,
               rho = rho,
               m = m,
               next_migrant = next_migrant)
  
  # run efficient C++ code
  output_raw <- SMC_MCMC_cpp(args)
  
  # TODO
  return(output_raw)
}

