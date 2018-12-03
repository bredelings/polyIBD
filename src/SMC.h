
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