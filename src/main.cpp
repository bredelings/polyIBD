
#include <Rcpp.h>
#include <RcppParallel.h>
#include <chrono>
#include "main.h"
#include "misc.h"
#include "probability.h"
#include "forwardBackward.h"

using namespace std;

//------------------------------------------------
// [[Rcpp::export]]
Rcpp::List runMCMC2_cpp(Rcpp::List args) {
    
    // start timer
    chrono::high_resolution_clock::time_point t1 = chrono::high_resolution_clock::now();
    
    // define parameters from inputs
    vector<int> x = Rcpp_to_vector_int(args["x"]);
    vector<double> p = Rcpp_to_vector_double(args["p"]);
    vector<int> SNP_dist = Rcpp_to_vector_int(args["SNP_dist"]);
    double rho = Rcpp_to_double(args["rho"]);
    double e1 = Rcpp_to_double(args["e1"]);
    double e2 = Rcpp_to_double(args["e2"]);
    int m_max = Rcpp_to_int(args["m_max"]);
    int burnin = Rcpp_to_int(args["burnin"]);
    int samples = Rcpp_to_int(args["samples"]);
    int L = x.size();
    
    // create lookup tables
    vector< vector< vector< vector< vector<double> > > > > emmission_lookup = define_emmission_lookup(m_max, L, p, e1, e2);
    vector< vector< vector<double> > > transition_lookup(L-1, vector< vector<double> >(m_max+1, vector<double> (m_max+1)));
    
    // initialise MCMC parameters
    int m1 = 1;
    int m2 = 1;
    int z_max = (m1>m2) ? m2 : m1;
    double f = 0.4;
    double logLike_old = 0;
    vector< vector<double> > frwrd_mat(m_max+1, vector<double>(L));
    vector< vector<double> > bkwrd_mat(m_max+1, vector<double>(L));
    vector< vector<double> > IBD_mat(m_max+1, vector<double>(L));
    
    // update transition probabilities
    update_transition_lookup(transition_lookup, f, rho, z_max);
    
    // call forward and backward algorithms
    logLike_old = forward_alg(frwrd_mat, x, emmission_lookup[m1-1][m2-1], transition_lookup, f);
    backward_alg(bkwrd_mat, x, emmission_lookup[m1-1][m2-1], transition_lookup);
    get_IBD(IBD_mat, frwrd_mat, bkwrd_mat);
    
    printMatrix(frwrd_mat);
    printMatrix(bkwrd_mat);
    printMatrix(IBD_mat);
    
    // end timer
    chrono::high_resolution_clock::time_point t2 = chrono::high_resolution_clock::now();
    chrono::duration<double> time_span = chrono::duration_cast< chrono::duration<double> >(t2-t1);
    print("   function completed in", time_span.count(), "seconds");
    
    // return values back to R
    return Rcpp::List::create(Rcpp::Named("IBD_mat")=IBD_mat);
}

