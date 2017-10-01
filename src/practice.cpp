
#include <Rcpp.h>
#include <RcppParallel.h>
#include "misc.h"
#include "probability.h"

using namespace std;

// [[Rcpp::export]]
Rcpp::List MCMC_VarMean(Rcpp::List args) {

    print("C code working");

    // NOTES
    // At the moment all our input arguments are in the form of an Rcpp::List. This is one of a special class of Rcpp objects, others include vectors and matrices. Although we could work directly with these Rcpp objects, I like to convert input arguments to standard C++ objects first (I'm just more familiar with how to use these). The overhead of converting objects should usually be small, but if not we might want to revisit this.

    // extract model parameters
    int n = Rcpp_to_int(args["n"]);
    double mu = Rcpp_to_double(args["mu"]);
    double sigma = Rcpp_to_double(args["sigma"]);

    // at this point we have all inputs as standard C++ objects, and we can proceed with the main logic of the code.

    // I like to keep track of execution time
    chrono::high_resolution_clock::time_point t1 = chrono::high_resolution_clock::now();

    // in this example we will create a vector of normal random variables.
    vector<double> x(n);
    for (int i=0; i<n; i++) {
        x[i] = rnorm1(mu,sigma);
    }

    // end timer
    chrono::high_resolution_clock::time_point t2 = chrono::high_resolution_clock::now();
    chrono::duration<double> time_span = chrono::duration_cast< chrono::duration<double> >(t2-t1);
    print("   function completed in", time_span.count(), "seconds");

    // return values back to R
    return Rcpp::List::create(Rcpp::Named("x")=x);
}
