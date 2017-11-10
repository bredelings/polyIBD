
#include <Rcpp.h>
#include <RcppParallel.h>
#include <chrono>
#include "main.h"
#include "misc.h"
#include "probability.h"

using namespace std;

//------------------------------------------------
// [[Rcpp::export]]
Rcpp::List runMCMC2_cpp(Rcpp::List args) {
    
    // start timer
    chrono::high_resolution_clock::time_point t1 = chrono::high_resolution_clock::now();
    
    // extract input arguments
    vector<int> x = Rcpp_to_vector_int(args["x"]);
    vector<double> p = Rcpp_to_vector_double(args["p"]);
    vector<int> SNP_dist = Rcpp_to_vector_int(args["SNP_dist"]);
    double rho = Rcpp_to_double(args["rho"]);
    int m_max = Rcpp_to_int(args["m_max"]);
    int burnin = Rcpp_to_int(args["burnin"]);
    int samples = Rcpp_to_int(args["samples"]);
    
    // define other parameters
    int L = x.size();
    
    // create lookup table
    vector< vector< vector< vector< vector<double> > > > > emmission_lookup = define_emmission_lookup(m_max, L, p);
    
    
    // end timer
    chrono::high_resolution_clock::time_point t2 = chrono::high_resolution_clock::now();
    chrono::duration<double> time_span = chrono::duration_cast< chrono::duration<double> >(t2-t1);
    print("   function completed in", time_span.count(), "seconds");
    
    // return values back to R
    return Rcpp::List::create(Rcpp::Named("x")=x);
}

//------------------------------------------------
// define emmission probability lookup table. This table is essentially a list over m1, then m2, then k, where k is the number of HMM states at each locus (i.e. k = 1 + the maximum level of IBD between samples). However, m2 must always be >= m1 in this list to minimise repetition. Therefore, actual values of m1 and m2 may need to be swapped when indexing this lookup table. Note, this also means that k is always equal to m1+1. At the final level of the list is a matrix with L rows and 9 columns giving emmission probabilities for each locus.
vector< vector< vector< vector< vector<double> > > > > define_emmission_lookup(int m_max, int L, vector<double> &p_vec) {
    
    // temporary vector
    vector<double> x(9);
    
    // loop through m1, m2 and k
    vector< vector< vector< vector< vector<double> > > > > emmission_lookup(m_max);
    for (int m1=1; m1<=m_max; m1++) {
        int m1_i = m1 - 1;
        emmission_lookup[m1_i] = vector< vector< vector< vector<double> > > >(m_max-m1_i);
        for (int m2=m1; m2<=m_max; m2++) {
            int m2_i = m2 - m1;
            int k = m1 + 1;
            emmission_lookup[m1_i][m2_i] = vector< vector< vector<double> > >(k);
            for (int z=0; z<k; z++) {
                emmission_lookup[m1_i][m2_i][z] = vector< vector<double> >(L, vector<double>(9));
                
                // loop through loci and HMM states
                for (int i=0; i<L; i++) {
                    double p = p_vec[i];
                    double q = 1-p;
                    
                    // calculate raw emmission probability
                    // if no IBD
                    if (z==0) {
                        
                        x[0] = pow(p,m1+m2);
                        x[1] = pow(p,m1) * (1 - pow(p,m2) - pow(q,m2));
                        x[2] = pow(p,m1) * pow(q,m2);
                        
                        x[3] = (1 - pow(p,m1) - pow(q,m1)) * pow(p,m2);
                        x[4] = (1 - pow(p,m1) - pow(q,m1)) * (1 - pow(p,m2) - pow(q,m2));
                        x[5] = (1 - pow(p,m1) - pow(q,m1)) * pow(q,m2);
                        
                        x[6] = pow(q,m1) * pow(p,m2);
                        x[7] = pow(q,m1) * (1 - pow(p,m2) - pow(q,m2));
                        x[8] = pow(q,m1+m2);
                        
                    }
                    // if some IBD
                    else {
                        
                        x[0] = pow(p,m1+m2-z);
                        x[1] = pow(p,m1) * (1 - pow(p,m2-z));
                        x[2] = 0;
                        
                        x[3] = pow(p,m2) * (1 - pow(p,m1-z));
                        x[4] = (1 - pow(p,m1) - pow(q,m1)) - pow(p,m2)*(1 - pow(p,m1-z)) - pow(q,m2)*(1 - pow(q,m1-z));
                        x[5] = pow(q,m2) * (1 - pow(q,m1-z));
                        
                        x[6] = 0;
                        x[7] = pow(q,m1) * (1 - pow(q,m2-z));
                        x[8] = pow(q,m1+m2-z);
                        
                    }
                    
                    emmission_lookup[m1_i][m2_i][z][i] = x;
                }
                
            }
        }
    }
    
    return(emmission_lookup);
}

