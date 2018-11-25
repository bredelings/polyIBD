
#include <Rcpp.h>
#include <RcppParallel.h>
#include <chrono>
#include "misc.h"
#include "probability.h"
#include "MCMC.h"
#include "SMC.h"

using namespace std;

//------------------------------------------------
// [[Rcpp::export]]
Rcpp::List runMCMC_cpp(Rcpp::List args, Rcpp::List args_functions) {
  
  // start timer
  chrono::high_resolution_clock::time_point t1 = chrono::high_resolution_clock::now();
  
  // define MCMC object
  MCMC mainMCMC(args, args_functions);
  
  // burn-in phase
  mainMCMC.burnin_MCMC(args_functions);
  
  // sampling phase
  mainMCMC.samp_MCMC(args_functions);
  
  // end timer
  chrono::high_resolution_clock::time_point t2 = chrono::high_resolution_clock::now();
  chrono::duration<double> time_span = chrono::duration_cast< chrono::duration<double> >(t2-t1);
  print("MCMC completed in", time_span.count(), "seconds");
  
  // return values back to R
  Rcpp::List ret;
  ret.push_back(Rcpp::wrap( mainMCMC.logLike_burnin_store ));
  ret.push_back(Rcpp::wrap( mainMCMC.logLike_store ));
  ret.push_back(Rcpp::wrap( mainMCMC.m1_store ));
  ret.push_back(Rcpp::wrap( mainMCMC.m2_store ));
  ret.push_back(Rcpp::wrap( mainMCMC.f_store ));
  ret.push_back(Rcpp::wrap( mainMCMC.f_ind_store ));
  ret.push_back(Rcpp::wrap( mainMCMC.k_store ));
  //ret.push_back(Rcpp::wrap( mainMCMC.sim_trans_n_store ));
  ret.push_back(Rcpp::wrap( mainMCMC.IBD_marginal ));
  ret.push_back(Rcpp::wrap( mainMCMC.accept_rate ));
  ret.push_back(Rcpp::wrap( time_span.count() ));
  
  Rcpp::StringVector ret_names;
  ret_names.push_back("logLike_burnin");
  ret_names.push_back("logLike");
  ret_names.push_back("m1");
  ret_names.push_back("m2");
  ret_names.push_back("f");
  ret_names.push_back("f_ind");
  ret_names.push_back("k");
//ret_names.push_back("sim_trans_n");
  ret_names.push_back("IBD_marginal");
  ret_names.push_back("accept_rate");
  ret_names.push_back("runTime");
  
  ret.names() = ret_names;
  return ret;
}

