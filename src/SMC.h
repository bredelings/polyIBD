
#pragma once

//------------------------------------------------
// draw from complete ancestral recombination graph under wright-fisher model
Rcpp::List wrightfisher_ARG_cpp(Rcpp::List args);

//------------------------------------------------
// draw from sequentially Markov coalescent under wright-fisher model
Rcpp::List wrightfisher_SMC_cpp(Rcpp::List args);

//------------------------------------------------
// draw from SMC under wright-fisher model by looping through discrete
// generations
Rcpp::List wrightfisher_SMC_naive_cpp(Rcpp::List args);

//------------------------------------------------
// draw from SMC under wright-fisher model by layering haplotypes one at a time
Rcpp::List wrightfisher_SMC_layer_cpp(Rcpp::List args);

//------------------------------------------------
// draw from complete ancestral recombination graph under continent-island model
Rcpp::List continentisland_ARG_cpp(Rcpp::List args);

//------------------------------------------------
// draw from SMC under continent-island model
Rcpp::List continentisland_SMC_cpp(Rcpp::List args);