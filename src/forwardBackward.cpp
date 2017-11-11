
#include <Rcpp.h>
#include "forwardBackward.h"
#include "misc.h"
#include "probability.h"

using namespace std;

//------------------------------------------------
// define emmission probability lookup table. This table is essentially a list over m1, then m2, then k, where k is the number of HMM states at each locus (i.e. k = 1 + the maximum level of IBD between samples). However, m2 must always be >= m1 in this list to minimise repetition. Therefore, actual values of m1 and m2 may need to be swapped when indexing this lookup table. Note, this also means that k is always equal to m1+1. At the final level of the list is a matrix with L rows and 9 columns giving emmission probabilities for each locus.
vector< vector< vector< vector< vector<double> > > > > define_emmission_lookup(int m_max, int L, vector<double> &p_vec, const double e1, const double e2) {
    
    // temporary vector for storing emmission probs
    vector<double> x_raw(9);
    vector<double> x_err(9);
    
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
                    
                    // calculate raw emmission probability (without error)
                    // if no IBD
                    if (z==0) {
                        
                        x_raw[0] = pow(p,m1+m2);
                        x_raw[1] = pow(p,m1) * (1 - pow(p,m2) - pow(q,m2));
                        x_raw[2] = pow(p,m1) * pow(q,m2);
                        
                        x_raw[3] = (1 - pow(p,m1) - pow(q,m1)) * pow(p,m2);
                        x_raw[4] = (1 - pow(p,m1) - pow(q,m1)) * (1 - pow(p,m2) - pow(q,m2));
                        x_raw[5] = (1 - pow(p,m1) - pow(q,m1)) * pow(q,m2);
                        
                        x_raw[6] = pow(q,m1) * pow(p,m2);
                        x_raw[7] = pow(q,m1) * (1 - pow(p,m2) - pow(q,m2));
                        x_raw[8] = pow(q,m1+m2);
                        
                    }
                    // if some IBD
                    else {
                        
                        x_raw[0] = pow(p,m1+m2-z);
                        x_raw[1] = pow(p,m1) * (1 - pow(p,m2-z));
                        x_raw[2] = 0;
                        
                        x_raw[3] = (1 - pow(p,m1-z)) * pow(p,m2);
                        x_raw[4] = (1 - pow(p,m1) - pow(q,m1)) - (1 - pow(p,m1-z))*pow(p,m2) - (1 - pow(q,m1-z))*pow(q,m2);
                        x_raw[5] = (1 - pow(q,m1-z)) * pow(q,m2);
                        
                        x_raw[6] = 0;
                        x_raw[7] = pow(q,m1) * (1 - pow(q,m2-z));
                        x_raw[8] = pow(q,m1+m2-z);
                        
                    }
                    
                    // incorporate error
                    x_err[0] = (1-e1)*((1-e1)*x_raw[0] + 0.5*e2*x_raw[1]) + 0.5*e2*((1-e1)*x_raw[3] + 0.5*e2*x_raw[4]);
                    x_err[1] = (1-e1)*((1-e2)*x_raw[1] + e1*(x_raw[0]+x_raw[2])) + 0.5*e2*((1-e2)*x_raw[4] + e1*(x_raw[3]+x_raw[5]));
                    x_err[2] = (1-e1)*((1-e1)*x_raw[2] + 0.5*e2*x_raw[1]) + 0.5*e2*((1-e1)*x_raw[5] + 0.5*e2*x_raw[4]);
                    
                    x_err[3] = (1-e2)*((1-e1)*x_raw[3] + 0.5*e2*x_raw[4]) + e1*((1-e1)*x_raw[0] + 0.5*e2*x_raw[1]) + e1*((1-e1)*x_raw[6] + 0.5*e2*x_raw[7]);
                    x_err[4] = (1-e2)*((1-e2)*x_raw[4] + e1*(x_raw[3]+x_raw[5])) + e1*((1-e2)*x_raw[1] + e1*(x_raw[0]+x_raw[2])) + e1*((1-e2)*x_raw[7] + e1*(x_raw[6]+x_raw[8]));
                    x_err[5] = (1-e2)*((1-e1)*x_raw[5] + 0.5*e2*x_raw[4]) + e1*((1-e1)*x_raw[2] + 0.5*e2*x_raw[1]) + e1*((1-e1)*x_raw[8] + 0.5*e2*x_raw[7]);
                    
                    x_err[6] = (1-e1)*((1-e1)*x_raw[6] + 0.5*e2*x_raw[7]) + 0.5*e2*((1-e1)*x_raw[3] + 0.5*e2*x_raw[4]);
                    x_err[7] = (1-e1)*((1-e2)*x_raw[7] + e1*(x_raw[6]+x_raw[8])) + 0.5*e2*((1-e2)*x_raw[4] + e1*(x_raw[3]+x_raw[5]));
                    x_err[8] = (1-e1)*((1-e1)*x_raw[8] + 0.5*e2*x_raw[7]) + 0.5*e2*((1-e1)*x_raw[5] + 0.5*e2*x_raw[4]);
                    
                    emmission_lookup[m1_i][m2_i][z][i] = x_err;
                }
                
            }
        }
    }
    
    return(emmission_lookup);
}

//------------------------------------------------
// update transition probability lookup table
void update_transition_lookup(vector< vector< vector<double> > > &transition_lookup, double f, double rho, int z_max) {
    
    // get dimensions
    int L = transition_lookup.size()+1;
    int s = transition_lookup[0].size();
    
    // define alpha from f
    double alpha = rho*f/(1-f);
    alpha = (alpha>OVERFLO) ? OVERFLO : alpha;
    
    // populate transition matrix
    for (int j=0; j<(L-1); j++) {
        for (int z1=0; z1<s; z1++) {
            fill(transition_lookup[j][z1].begin(), transition_lookup[j][z1].end(), 0);
        }
        for (int z1=0; z1<=z_max; z1++) {
            for (int z2=0; z2<=z_max; z2++) {
                if (z1==z2) {
                    transition_lookup[j][z1][z2] = 1000/double(1000 + z_max);
                } else {
                    transition_lookup[j][z1][z2] = 1/double(1000 + z_max);
                }
            }
        }
    }
    
}

//------------------------------------------------
// text
double forward_alg(vector< vector<double> > &frwrd_mat, const vector<int> &x, const vector< vector< vector<double> > > &emmission_lookup, const vector< vector< vector<double> > > &transition_lookup, double f) {
    
    // get basic parameters from inputs
    int L = x.size();
    int z_max = emmission_lookup.size()-1;
    print(z_max);
    
    // clear frwrd_mat
    for (int i=0; i<int(frwrd_mat.size()); i++) {
        fill(frwrd_mat[i].begin(), frwrd_mat[i].end(), 0);
    }
    
    // carry out first step of algorithm
    double frwrd_sum = 0;
    double logLike = 0;
    for (int z=0; z<=z_max; z++) {
        frwrd_mat[z][0] = R::dbinom(z,z_max,f,false) * emmission_lookup[z][0][x[0]];
        frwrd_sum += frwrd_mat[z][0];
    }
    logLike += log(frwrd_sum);
    for (int z=0; z<=z_max; z++) {
        frwrd_mat[z][0] /= frwrd_sum;
    }
    
    // carry out remaining steps of algorithm
    for (int j=1; j<L; j++) {
        frwrd_sum = 0;
        for (int z=0; z<=z_max; z++) {
            for (int i=0; i<=z_max; i++) {
                frwrd_mat[z][j] += frwrd_mat[i][j-1] * transition_lookup[j-1][i][z];
            }
            frwrd_mat[z][j] *= emmission_lookup[z][j][x[j]];
            frwrd_sum += frwrd_mat[z][j];
        }
        logLike += log(frwrd_sum);
        for (int z=0; z<=z_max; z++) {
            frwrd_mat[z][j] /= frwrd_sum;
        }
    }
    
    return(logLike);
}

//------------------------------------------------
// text
void backward_alg(vector< vector<double> > &bkwrd_mat, const vector<int> &x, const vector< vector< vector<double> > > &emmission_lookup, const vector< vector< vector<double> > > &transition_lookup) {
    
    // get basic parameters from inputs
    int L = x.size();
    int z_max = emmission_lookup.size()-1;
    
    // reset bkwrd_mat
    for (int i=0; i<int(bkwrd_mat.size()); i++) {
        fill(bkwrd_mat[i].begin(), bkwrd_mat[i].end(), 1);
    }
    
    // loop through loci backwards from penultimate
    double bkwrd_sum = 0;
    for (int j=(L-2); j>=0; j--) {
        bkwrd_sum = 0;
        for (int z=0; z<=z_max; z++) {
            for (int i=0; i<=z_max; i++) {
                bkwrd_mat[z][j] += transition_lookup[j][z][i]*bkwrd_mat[i][j+1]*emmission_lookup[i][j+1][x[j+1]];
            }
            bkwrd_sum += bkwrd_mat[z][j];
        }
        for (int z=0; z<=z_max; z++) {
            bkwrd_mat[z][j] /= bkwrd_sum;
        }
    }
    
}

//------------------------------------------------
// text
void get_IBD(vector< vector<double> > &IBD_mat, const vector< vector<double> > &frwrd_mat, const vector< vector<double> > &bkwrd_mat) {
    
    int L = IBD_mat[0].size();
    int s = IBD_mat.size();
    
    double IBD_sum = 0;
    for (int j=0; j<L; j++) {
        IBD_sum = 0;
        for (int z=0; z<s; z++) {
            IBD_mat[z][j] = frwrd_mat[z][j] * bkwrd_mat[z][j];
            IBD_sum += IBD_mat[z][j];
        }
        for (int z=0; z<s; z++) {
            IBD_mat[z][j] /= IBD_sum;
        }
    }
    
}


