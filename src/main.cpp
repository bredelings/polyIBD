
#include <Rcpp.h>
#include <RcppParallel.h>
#include <chrono>
#include "main.h"
#include "misc.h"
#include "probability.h"
#include "MCMC.h"

using namespace std;

//------------------------------------------------
// [[Rcpp::export]]
Rcpp::List runMCMC2_cpp(Rcpp::List args, Rcpp::List args_functions) {
    
    // start timer
    chrono::high_resolution_clock::time_point t1 = chrono::high_resolution_clock::now();
    
    // define MCMC object
    MCMC mainMCMC(args, args_functions);
    
    // burn-in phase
    mainMCMC.burnin_MCMC(args_functions);
    
    // sampling phase
    mainMCMC.run_MCMC(args_functions);
    
    // end timer
    chrono::high_resolution_clock::time_point t2 = chrono::high_resolution_clock::now();
    chrono::duration<double> time_span = chrono::duration_cast< chrono::duration<double> >(t2-t1);
    print("MCMC completed in", time_span.count(), "seconds");
    
    // return values back to R
    return Rcpp::List::create(Rcpp::Named("logLike_burnin")=mainMCMC.logLike_burnin_store, Rcpp::Named("logLike")=mainMCMC.logLike_store, Rcpp::Named("m1")=mainMCMC.m1_store, Rcpp::Named("m2")=mainMCMC.m2_store, Rcpp::Named("f")=mainMCMC.f_store, Rcpp::Named("rho")=mainMCMC.rho_store, Rcpp::Named("IBD_marginal")=mainMCMC.IBD_marginal);
}

