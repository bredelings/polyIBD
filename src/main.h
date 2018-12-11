
#include <Rcpp.h>

#pragma once

//------------------------------------------------
// draw from complete ancestral recombination graph under continent-island model
Rcpp::List continentisland_ARG_cpp(Rcpp::List args);

//------------------------------------------------
// draw from conditional SMC under continent-island model
Rcpp::List continentisland_SMC_conditional_cpp(Rcpp::List args);

//------------------------------------------------
// MCMC under continent-island model
Rcpp::List SMC_MCMC_cpp(Rcpp::List args);

