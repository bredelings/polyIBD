
#include <Rcpp.h>
#include <RcppParallel.h>
#include "misc.h"
#include "probability.h"

using namespace std;

// [[Rcpp::export]]
Rcpp::List MCMC_VarMean_cpp(Rcpp::List args) {

    print("C code working");



    // Start Time
    chrono::high_resolution_clock::time_point t1 = chrono::high_resolution_clock::now();

    // Input from R
    std::vector<double> inputFromR = Rcpp_to_vector_double(args["input"]);

    // Define Prior Distribution
    int priordistributioner (int mu, int sigma)
    {
    	int muprior
		int sigmaprior


    }






    // end timer
    chrono::high_resolution_clock::time_point t2 = chrono::high_resolution_clock::now();
    chrono::duration<double> time_span = chrono::duration_cast< chrono::duration<double> >(t2-t1);
    print("   function completed in", time_span.count(), "seconds");

    // return values back to R
    return Rcpp::List::create(Rcpp::Named("chain")=inputFromR);
}
