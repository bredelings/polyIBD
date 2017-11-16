
#include <Rcpp.h>
#include <RcppParallel.h>
#include "MCMC.h"
#include "misc.h"
#include "probability.h"

using namespace std;

//------------------------------------------------
// MCMC::
// constructor for MCMC class
MCMC::MCMC(Rcpp::List args, Rcpp::List args_functions) {
    
    // data and model parameters
    x = Rcpp_to_vector_int(args["x"]);
    L = x.size();
    p_vec = Rcpp_to_vector_double(args["p"]);
    SNP_dist = Rcpp_to_vector_int(args["SNP_dist"]);
    e1 = Rcpp_to_double(args["e1"]);
    e2 = Rcpp_to_double(args["e2"]);
    m_max = Rcpp_to_int(args["m_max"]);
    
    // MCMC parameters
    burnin = Rcpp_to_int(args["burnin"]);
    samples = Rcpp_to_int(args["samples"]);
    reportIteration = Rcpp_to_int(args["reportIteration"]);
    
    // define lookup tables
    define_emmission_lookup();
    transition_lookup = vector< vector< vector<double> > >(L-1, vector< vector<double> >(m_max+1, vector<double> (m_max+1)));
    
    // initialise transient MCMC objects
    m1 = 1;
    m2 = 1;
    z_max = (m1<m2) ? m1 : m2;
    f = 0.01;
    rho = 0.01;
    logLike_old = 0;
    frwrd_mat = vector< vector<double> >(m_max+1, vector<double>(L));
    bkwrd_mat = vector< vector<double> >(m_max+1, vector<double>(L));
    IBD_mat = vector< vector<double> >(m_max+1, vector<double>(L));
    
    // first pass through forward-backward algorithm
    update_transition_lookup(f, rho, m1, m2, args_functions["getTransProbs"]);
    logLike_old = forward_alg(m1, m2);
    backward_alg(m1, m2);
    get_IBD();
    
    // objects for storing MCMC results
    logLike_burnin_store = vector<double>(burnin);
    logLike_store = vector<double>(samples);
    m1_store = vector<int>(samples);
    m2_store = vector<int>(samples);
    f_store = vector<double>(samples);
    rho_store = vector<double>(samples);
    IBD_store = vector< vector< vector<double> > >(samples);
    IBD_weight = vector<double>(samples);
    IBD_marginal = vector< vector<double> >(m_max+1, vector<double>(L));
    
    logLike_burnin_store[0] = logLike_old;
    IBD_store[0] = IBD_mat;
    IBD_weight[0] = 1;
    
    // misc objects
    m1_weight_stay = 1;
    m1_weight_move = 1;
    m2_weight_stay = 1;
    m2_weight_move = 1;
    f_propSD = 0.1;
    rho_propSD = 0.1;
    IBD_index = 0;
    
}

//------------------------------------------------
// MCMC::
// run MCMC
void MCMC::burnin_MCMC(Rcpp::List args_functions) {
    
    print("Running MCMC");
    print("   burnin phase");
    for (int rep=1; rep<burnin; rep++) {
        
        // report progress
        if (reportIteration>0) {
            if (((rep+1) % reportIteration)==0) {
                print("      iteration",rep+1);
            }
        }
        
        // propose m1 and m2
        int m1_prop = propose_m(m1, m1_weight_move, m1_weight_stay);
        int m2_prop = propose_m(m2, m2_weight_move, m2_weight_stay);
        
        // if no change in m, propose either f or rho
        double f_prop = f;
        double rho_prop = rho;
        if (m1_prop==m1 && m2_prop==m2) {
            if (rbernoulli1(0.5)) {
                f_prop = rnorm1_interval(f, f_propSD, 0, 1);
            } else {
                rho_prop = rnorm1_interval(rho, rho_propSD, 0, 0.1);
            }
        }
        
        // update transition probabilities and calculate new likelihood
        update_transition_lookup(f_prop, rho_prop, m1_prop, m2_prop, args_functions["getTransProbs"]);
        double logLike_new = forward_alg(m1_prop, m2_prop);
        
        // Metropolis-Hastings step
        // if accept
        if (log(runif_0_1()) < (logLike_new-logLike_old)) {
            
            // update m1 and m2 sampling weights
            m1_weight_move = (m1==m1_prop) ? m1_weight_move : ++m1_weight_move;
            m2_weight_move = (m2==m2_prop) ? m2_weight_move : ++m2_weight_move;
            
            // update f_propSD
            if (m1==m1_prop && m2==m2_prop && f_prop!=f) {
                f_propSD  += (1-0.23)/sqrt(double(rep));
            }
            
            // update rho_propSD by Robbins-Monro
            if (m1==m1_prop && m2==m2_prop && rho_prop!=rho) {
                rho_propSD  += (1-0.23)/sqrt(double(rep));
            }
            
            // update parameter values and likelihood
            f = f_prop;
            rho = rho_prop;
            m1 = m1_prop;
            m2 = m2_prop;
            logLike_old = logLike_new;
            
        }
        // if reject
        else {
            
            // update m1 and m2 sampling weights
            m1_weight_stay = (m1==m1_prop) ? m1_weight_stay : ++m1_weight_stay;
            m2_weight_stay = (m2==m2_prop) ? m2_weight_stay : ++m2_weight_stay;
            
            // update f_propSD by Robbins-Monro
            if (m1==m1_prop && m2==m2_prop && f_prop!=f) {
                f_propSD  -= 0.23/sqrt(double(rep));
                f_propSD = (f_propSD < 0) ? -f_propSD : f_propSD;
            }
            
            // update rho_propSD by Robbins-Monro
            if (m1==m1_prop && m2==m2_prop && rho_prop!=rho) {
                rho_propSD  -= 0.23/sqrt(double(rep));
                rho_propSD = (rho_propSD < 0) ? -rho_propSD : rho_propSD;
            }
            
        }
        
        // store logLike
        logLike_burnin_store[rep] = logLike_old;
        
    }   // end MCMC loop
    
}

//------------------------------------------------
// MCMC::
// run MCMC
void MCMC::run_MCMC(Rcpp::List args_functions) {
    
    print("   sampling phase");
    IBD_index = 0;
    for (int rep=0; rep<samples; rep++) {
        
        // report progress
        if (reportIteration>0) {
            if (((rep+1) % reportIteration)==0) {
                print("      iteration",rep+1);
            }
        }
        
        // propose m1 and m2
        int m1_prop = propose_m(m1, m1_weight_move, m1_weight_stay);
        int m2_prop = propose_m(m2, m2_weight_move, m2_weight_stay);
        
        // if no change in m, propose either f or rho
        double f_prop = f;
        double rho_prop = rho;
        if (m1_prop==m1 && m2_prop==m2) {
            if (rbernoulli1(0.5)) {
                f_prop = rnorm1_interval(f, f_propSD, 0, 1);
            } else {
                rho_prop = rnorm1_interval(rho, rho_propSD, 0, 0.1);
            }
        }
        
        // update transition probabilities and calculate new likelihood
        update_transition_lookup(f_prop, rho_prop, m1_prop, m2_prop, args_functions["getTransProbs"]);
        double logLike_new = forward_alg(m1_prop, m2_prop);
        
        // Metropolis-Hastings step
        if (log(runif_0_1()) < (logLike_new-logLike_old)) {
            
            // update parameter values and likelihood
            f = f_prop;
            rho = rho_prop;
            m1 = m1_prop;
            m2 = m2_prop;
            logLike_old = logLike_new;
            
            // update and store IBD_mat
            backward_alg(m1, m2);
            get_IBD();
            IBD_index++;
            IBD_store[IBD_index] = IBD_mat;
            
        }
        
        // increase current IBD_weight
        IBD_weight[IBD_index] ++;
        
        // store current values
        logLike_store[rep] = logLike_old;
        m1_store[rep] = m1;
        m2_store[rep] = m2;
        f_store[rep] = f;
        rho_store[rep] = rho;
        
    }   // end MCMC loop
    
    
    // calculate marginal IBD matrix over all MCMC output
    for (int rep=0; rep<samples; rep++) {
        if (IBD_weight[rep]==0) {
            break;
        }
        for (int i=0; i<(m_max+1); i++) {
            for (int j=0; j<L; j++) {
                IBD_marginal[i][j] += IBD_weight[rep]*IBD_store[rep][i][j];
            }
        }
    }
    
}

//------------------------------------------------
// MCMC::
// define emmission probability lookup table. This table is essentially a list over m1, then m2, then k, where k is the number of HMM states at each locus (i.e. k = 1 + the maximum level of IBD between samples). At the final level of the list is a matrix with L rows and 9 columns giving emmission probabilities for each locus.
void MCMC::define_emmission_lookup() {
    
    // temporary vector for storing emmission probs
    vector<double> x_raw(9);
    vector<double> x_err(9);
    
    // loop through m1, m2 and k
    emmission_lookup = vector< vector< vector< vector< vector<double> > > > >(m_max);
    for (int m1=1; m1<=m_max; m1++) {
        int m1_i = m1 - 1;
        emmission_lookup[m1_i] = vector< vector< vector< vector<double> > > >(m_max);
        for (int m2=1; m2<=m_max; m2++) {
            int m2_i = m2 - 1;
            int k = m1 + 1;
            k = (m1<m2) ? k : m2 + 1;
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
}

//------------------------------------------------
// MCMC::
// update transition probability lookup table
void MCMC::update_transition_lookup(double f, double rho, int m1, int m2, Rcpp::Function getTransProbs) {
    
    // get z_max
    int z_max = (m1 < m2) ? m1 : m2;
    
    // get eigenvalue solutions to rate matrix
    Rcpp::List Elist = getTransProbs(f, rho, z_max);
    vector<double> Evalues = Rcpp_to_vector_double(Elist["Evalues"]);
    vector< vector<double> > Evectors = Rcpp_to_mat_double(Elist["Evectors"]);
    vector< vector<double> > Esolve = Rcpp_to_mat_double(Elist["Esolve"]);
    
    // populate transition matrix
    for (int j=0; j<(L-1); j++) {
        // clear existing values
        for (int z1=0; z1<(m_max+1); z1++) {
            fill(transition_lookup[j][z1].begin(), transition_lookup[j][z1].end(), 0);
        }
        // populate based on rate matrix solution
        for (int z1=0; z1<(z_max+1); z1++) {
            for (int z2=0; z2<(z_max+1); z2++) {
                for (int i=0; i<(z_max+1); i++) {
                    transition_lookup[j][z1][z2] += Evectors[z2][i]*Esolve[i][z1] * exp(Evalues[i] * SNP_dist[j]);
                }
            }
        }
    }
    
}

//------------------------------------------------
// MCMC::
// text
double MCMC::forward_alg(int m1, int m2) {
    
    // get z_max
    int z_max = (m1 < m2) ? m1 : m2;
    
    // clear frwrd_mat
    for (int i=0; i<int(frwrd_mat.size()); i++) {
        fill(frwrd_mat[i].begin(), frwrd_mat[i].end(), 0);
    }
    
    // carry out first step of algorithm
    double frwrd_sum = 0;
    double logLike = 0;
    for (int z=0; z<(z_max+1); z++) {
        frwrd_mat[z][0] = R::dbinom(z,z_max,f,false) * emmission_lookup[m1-1][m2-1][z][0][x[0]];
        frwrd_sum += frwrd_mat[z][0];
    }
    logLike += log(frwrd_sum);
    for (int z=0; z<(z_max+1); z++) {
        frwrd_mat[z][0] /= frwrd_sum;
    }
    
    // carry out remaining steps of algorithm
    for (int j=1; j<L; j++) {
        frwrd_sum = 0;
        for (int z=0; z<(z_max+1); z++) {
            for (int i=0; i<(z_max+1); i++) {
                frwrd_mat[z][j] += frwrd_mat[i][j-1] * transition_lookup[j-1][i][z];
            }
            frwrd_mat[z][j] *= emmission_lookup[m1-1][m2-1][z][j][x[j]];
            frwrd_sum += frwrd_mat[z][j];
        }
        logLike += log(frwrd_sum);
        for (int z=0; z<(z_max+1); z++) {
            frwrd_mat[z][j] /= frwrd_sum;
        }
    }
    
    return(logLike);
}

//------------------------------------------------
// MCMC::
// text
void MCMC::backward_alg(int m1, int m2) {
    
    // get z_max
    int z_max = (m1 < m2) ? m1 : m2;
    
    // clear bkwrd_mat, set final value to 1
    for (int i=0; i<int(bkwrd_mat.size()); i++) {
        fill(bkwrd_mat[i].begin(), bkwrd_mat[i].end(), 0);
        bkwrd_mat[i][L-1] = 1;
    }
    
    // loop through loci backwards from penultimate
    double bkwrd_sum = 0;
    for (int j=(L-2); j>=0; j--) {
        bkwrd_sum = 0;
        for (int z=0; z<(z_max+1); z++) {
            for (int i=0; i<(z_max+1); i++) {
                bkwrd_mat[z][j] += transition_lookup[j][z][i]*bkwrd_mat[i][j+1]*emmission_lookup[m1-1][m2-1][i][j+1][x[j+1]];
            }
            bkwrd_sum += bkwrd_mat[z][j];
        }
        for (int z=0; z<(z_max+1); z++) {
            bkwrd_mat[z][j] /= bkwrd_sum;
        }
    }
    
}

//------------------------------------------------
// MCMC::
// text
void MCMC::get_IBD() {
    
    // get z_max
    int z_max = (m1 < m2) ? m1 : m2;
    
    // take product of forward and backward matrices, and normalise
    double IBD_sum = 0;
    for (int j=0; j<L; j++) {
        IBD_sum = 0;
        for (int z=0; z<(z_max+1); z++) {
            IBD_mat[z][j] = frwrd_mat[z][j] * bkwrd_mat[z][j];
            IBD_sum += IBD_mat[z][j];
        }
        for (int z=0; z<(z_max+1); z++) {
            IBD_mat[z][j] /= IBD_sum;
        }
    }
    
}

//------------------------------------------------
// MCMC::
// text
double MCMC::propose_m(double m_current, double weight_move, double weight_stay) {
    
    double m_prop = m_current;
    if (rbernoulli1( weight_move/double(weight_move + weight_stay) )) {
        m_prop = rbernoulli1(0.5) ? m_current+1 : m_current-1;
        m_prop = (m_prop==0 || m_prop>m_max) ? m_current : m_prop;
    }
    return(m_prop);
    
}

