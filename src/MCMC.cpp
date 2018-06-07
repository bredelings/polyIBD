
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
  rho = Rcpp_to_double(args["rho"]);
  k_max = Rcpp_to_double(args["k_max"]);
  SNP_dist = Rcpp_to_vector_int(args["SNP_dist"]);
  e1 = Rcpp_to_double(args["e1"]);
  e2 = Rcpp_to_double(args["e2"]);
  m_max = Rcpp_to_int(args["m_max"]);
  
  // MCMC parameters
  burnin = Rcpp_to_int(args["burnin"]);
  samples = Rcpp_to_int(args["samples"]);
  reportIteration = Rcpp_to_int(args["reportIteration"]);
  
  // define lookup tables
  // The emmission lookup table is fully defined here. The transition lookup table is defined as empty, and will be updated with new values throughout the MCMC.
  define_emmission_lookup();
  transition_lookup = vector< vector< vector<double> > >(L-1, vector< vector<double> >(m_max+1, vector<double> (m_max+1)));
  
  // initialise transient MCMC objects
  m1 = 1;
  m2 = 1;
  z_max = (m1<m2) ? m1 : m2;
  f = 0.01;
  k = 1;
  k_max = 50;
  logLike_old = 0;
  frwrd_mat = vector< vector< double> >(m_max+1, vector< double>(L));
  bkwrd_mat = vector< vector< double> >(m_max+1, vector< double>(L));
  IBD_mat = vector< vector< double> >(m_max+1, vector<double>(L));
  
  // objects for storing MCMC results
  logLike_burnin_store = vector< double>(burnin);
  logLike_store = vector< double>(samples);
  m1_store = vector<int>(samples);
  m2_store = vector<int>(samples);
  f_store = vector<double>(samples);
  k_store = vector<double>(samples);
  IBD_marginal = vector< vector< double> >(m_max+1, vector< double>(L));
  accept_rate = 0;
  
  // temp objects
  f_ind = 0;
  f_ind_store = vector<double>(samples);
  // sim_trans_n = 0;
  // sim_trans_n_store = vector<int>(samples);
  
  // calculate initial likelihood
  update_transition_lookup(f, rho, k, m1, m2, args_functions["getTransProbs"]);
  
  
  logLike_old = forward_alg(m1, m2);
  logLike_burnin_store[0] = logLike_old;

  // misc objects
  // m weights dictate the chance of proposing a new m value ("move") vs. sticking with the current value ("stay"). These weights are updated during the burn-in phase, converging to an efficient proposal distribution. The chance of staying is never allowed to exceed 100* the chance of moving, to ensure there is always some chance of exploring new m values.
  // f_propSD and k_propSD are updated during the burn-in phase using Robbins-Monro (although k_propSD isn't being used right now).
  m1_weight_stay = 1;
  m1_weight_move = 1;
  m2_weight_stay = 1;
  m2_weight_move = 1;
  f_propSD = 0.2;
  k_weight_stay = 1;
  k_weight_move = 1;
}

//------------------------------------------------
// MCMC::
// Run burn-in phase of MCMC. Each iteration we alternate which parameters are updated. First, we give m1 and m2 a chance to update, based on the sampling weights. If m1 or m2 do change then we don't update the other parameters. If by chance m1 and m2 stay the same we choose one of the other parameters to update with equal probability. The advantage to this method is that, because only one parameter is updated each iteration, we can only change the proposal standard deviation relating to that parameter (by Robbins-Monro). If multiple parameters changed each iteration we could not change these proposal standard deviations independently.
void MCMC::burnin_MCMC(Rcpp::List args_functions) {
  // announce burn-in phase
  print("Running MCMC");
  print("   burnin phase");
  
  // loop through burn-in iterations
  for (int rep=1; rep<burnin; rep++) {    // start in second rep
    
    // report progress
    if (reportIteration>0) {
      if (((rep+1) % reportIteration)==0) {
        print("      iteration",rep+1);
      }
    }
    
    // propose m1 and m2
    int m1_prop = propose_m(m1, m1_weight_move, m1_weight_stay);
    int m2_prop = propose_m(m2, m2_weight_move, m2_weight_stay);
    
    // if no change in m then propose either f or k
    double f_prop = f;
    int k_prop = k;
    if (m1_prop==m1 && m2_prop==m2) {
      if (rbernoulli1(0.5)) {
        f_prop = rnorm1_interval(f, f_propSD, 0, 1);
        
      } else {
        //some proposal distribution for k -- consider weights later
        k_prop = propose_k(k, k_weight_move, k_weight_stay);
      }
    }
    
    // update transition probabilities and calculate new likelihood
    update_transition_lookup(f_prop, rho, k_prop, m1_prop, m2_prop, args_functions["getTransProbs"]);
    double logLike_new = forward_alg(m1_prop, m2_prop);
    
    
    
    // Metropolis-Hastings step
    // note that all proposal distributions are symmetric, therefore no Hastings step is required
    // if accept
    if (log(runif_0_1()) < (logLike_new-logLike_old)) {
      
      // update m1 and m2 sampling weights
      m1_weight_move = (m1==m1_prop) ? m1_weight_move : ++m1_weight_move;
      m2_weight_move = (m2==m2_prop) ? m2_weight_move : ++m2_weight_move;
      
      // or update f_propSD
      if (m1==m1_prop && m2==m2_prop && f_prop!=f) {
        f_propSD  += (1-0.23)/sqrt(double(rep));
      }
      
      // or update k sample weights
      if (m1==m1_prop && m2==m2_prop && k_prop!=k) {
        k_weight_move = (k==k_prop) ? k_weight_move : ++k_weight_move;
      }
      
      // update parameter values and likelihood
      f = f_prop;
      k = k_prop;
      m1 = m1_prop;
      m2 = m2_prop;
      logLike_old = logLike_new;
      
    }
    // if reject
    else {
      
      // update m1 and m2 sampling weights
      m1_weight_stay = (m1==m1_prop) ? m1_weight_stay : ++m1_weight_stay;
      m2_weight_stay = (m2==m2_prop) ? m2_weight_stay : ++m2_weight_stay;
      
      // limit m1 and m2 sampling weights
      m1_weight_stay = (m1_weight_stay > 100*m1_weight_move) ? 100*m1_weight_move : m1_weight_stay;
      m2_weight_stay = (m2_weight_stay > 100*m2_weight_move) ? 100*m2_weight_move : m2_weight_stay;
      
      // update f_propSD
      if (m1==m1_prop && m2==m2_prop && f_prop!=f) {
        f_propSD  -= 0.23/sqrt(double(rep));
        f_propSD = (f_propSD < 0) ? -f_propSD : f_propSD;
      }
      
      // update k standard deviation
      if (m1==m1_prop && m2==m2_prop && k_prop!=k) {
        k_weight_stay = (k==k_prop) ? k_weight_stay : ++k_weight_stay;
        // limit k sampling weights
        k_weight_stay = (k_weight_stay > 100*k_weight_move) ? 100*k_weight_move : k_weight_stay;
      }
      
    }
    // store logLike
    logLike_burnin_store[rep] = logLike_old;
    
  }   // end MCMC loop
  
}

//------------------------------------------------
// MCMC::
// Run sampling phase of MCMC. Most of this function is very similar to the burn-in phase. Differences are that sampling weights and standard deviations are not updated by Robbins-Monro, and a greater number of outputs are stored.
void MCMC::run_MCMC(Rcpp::List args_functions) {
  
  // announce burn-in phase
  print("   sampling phase");
  
  // calculate initial IBD matrix
  backward_alg(m1, m2);
  get_IBD();
  
  // loop through sampling iterations
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
    
    // if no change in m, therefore propose either f or k
    double f_prop = f;
    int k_prop = k;
    if (m1_prop==m1 && m2_prop==m2) {
      if (rbernoulli1(0.5)) {
        f_prop = rnorm1_interval(f, f_propSD, 0, 1);
      } else {
        k_prop = propose_k(k, k_weight_move, k_weight_stay);
      }
    }
    
    // update transition probabilities and calculate new likelihood
    update_transition_lookup(f_prop, rho, k_prop, m1_prop, m2_prop, args_functions["getTransProbs"]);
    double logLike_new = forward_alg(m1_prop, m2_prop);
    
    // Metropolis-Hastings step
    // note that all proposal distributions are symmetric, therefore no Hastings step is required
    // if accept
    if (log(runif_0_1()) < (logLike_new-logLike_old)) {
      
      // update parameter values and likelihood
      f = f_prop;
      k = k_prop;
      m1 = m1_prop;
      m2 = m2_prop;
      logLike_old = logLike_new;
      
      // update IBD_mat
      backward_alg(m1, m2);
      get_IBD();
      
      // update acceptance rate
      accept_rate ++;
    }
    
    // store current values
    logLike_store[rep] = logLike_old;
    m1_store[rep] = m1;
    m2_store[rep] = m2;
    f_store[rep] = f;
    f_ind_store[rep] = f_ind;
    k_store[rep] = k;
    //sim_trans_n_store[rep] = sim_trans_n;
    
    // add current IBD_mat to running estimate
    for (int i=0; i<(m_max+1); i++) {
      for (int j=0; j<L; j++) {
        IBD_marginal[i][j] += IBD_mat[i][j];
      }
    }
    
  }   // end MCMC loop
  
  // finalise IBD_mat
  for (int i=0; i<(m_max+1); i++) {
    for (int j=0; j<L; j++) {
      IBD_marginal[i][j] /= double(samples);
    }
  }
}

//------------------------------------------------
// MCMC::
// define emmission probability lookup table. This table is essentially a list over m1, then m2, then k, where k is the number of HMM states at each locus (i.e. k = 1 + the maximum level of IBD between samples). At the final level of the list is a matrix with L rows and 16 columns giving emmission probabilities for each locus. The 16 columns correspond to the 16 possible genotypes combinations (NA means missing data):
// 0 = {NA, NA}
// 1 = {NA, A}
// 2 = {NA, Aa}
// 3 = {NA, a}
// 4 = {A, NA}
// 5 = {A,A}
// 6 = {A,Aa}
// 7 = {A,a}
// 8 = {Aa, NA}
// 9 = {Aa,A}
// 10 = {Aa,Aa}
// 11 = {Aa,a}
// 12 = {a, NA}
// 13 = {a,A}
// 14 = {a,Aa}
// 15 = {a,a}

void MCMC::define_emmission_lookup() {
  
  // temporary vector for storing emmission probs
  vector<double> x_raw(16);
  vector<double> x_err(16);
  
  // loop through m1
  emmission_lookup = vector< vector< vector< vector< vector<double> > > > >(m_max);
  for (int m1=1; m1<=m_max; m1++) {
    
    // loop through m2
    emmission_lookup[m1-1] = vector< vector< vector< vector<double> > > >(m_max);
    for (int m2=1; m2<=m_max; m2++) {
      
      // define k as 1 + maximum IBD
      int k = m1 + 1;
      k = (m1<m2) ? k : m2 + 1;
      
      // loop through z
      emmission_lookup[m1-1][m2-1] = vector< vector< vector<double> > >(k);
      for (int z=0; z<k; z++) {
        emmission_lookup[m1-1][m2-1][z] = vector< vector<double> >(L, vector<double>(16));
        
        // loop through loci
        for (int i=0; i<L; i++) {
          double p = p_vec[i];
          double q = 1-p;
          
          // calculate raw emmission probability (without error)
          // if no IBD
          if (z==0) {
            
            x_raw[0] = 1;
            x_raw[1] = pow(p,m2);
            x_raw[2] = (1 - pow(p,m2) - pow(q,m2));
            x_raw[3] = pow(q,m2);
            
            x_raw[4] = pow(p,m1);
            x_raw[5] = pow(p,m1+m2);
            x_raw[6] = pow(p,m1) * (1 - pow(p,m2) - pow(q,m2));
            x_raw[7] = pow(p,m1) * pow(q,m2);
            
            x_raw[8] = (1 - pow(p,m1) - pow(q,m1));
            x_raw[9] = (1 - pow(p,m1) - pow(q,m1)) * pow(p,m2);
            x_raw[10] = (1 - pow(p,m1) - pow(q,m1)) * (1 - pow(p,m2) - pow(q,m2));
            x_raw[11] = (1 - pow(p,m1) - pow(q,m1)) * pow(q,m2);
            
            x_raw[12] = pow(q,m1);
            x_raw[13] = pow(q,m1) * pow(p,m2);
            x_raw[14] = pow(q,m1) * (1 - pow(p,m2) - pow(q,m2));
            x_raw[15] = pow(q,m1+m2);
            
          }
          // if some IBD
          else {
            
            x_raw[0] = 1;
            x_raw[1] = pow(p,m2);
            x_raw[2] = (1 - pow(p,m2) - pow(q,m2));
            x_raw[3] = pow(q,m2);
            
            x_raw[4] = pow(p,m1);
            x_raw[5] = pow(p,m1+m2-z);
            x_raw[6] = pow(p,m1) * (1 - pow(p,m2-z));
            x_raw[7] = 0;
            
            x_raw[8] = (1 - pow(p,m1) - pow(q,m1));
            x_raw[9] = (1 - pow(p,m1-z)) * pow(p,m2);
            x_raw[10] = (1 - pow(p,m1) - pow(q,m1)) - (1 - pow(p,m1-z))*pow(p,m2) - (1 - pow(q,m1-z))*pow(q,m2);
            x_raw[11] = (1 - pow(q,m1-z)) * pow(q,m2);
            
            x_raw[12] = pow(q,m1);
            x_raw[13] = 0;
            x_raw[14] = pow(q,m1) * (1 - pow(q,m2-z));
            x_raw[15] = pow(q,m1+m2-z);
            
          }
          
          // incorporate error
          // e1 = prob homo incorrectly called as het
          // e2 = prob het incorrectly called as homo
          x_err[0] = x_raw[0];
          x_err[1] = (1-e1)*x_raw[1] + 0.5*e2*x_raw[2];
          x_err[2] = (1-e2)*x_raw[2] + e1*(x_raw[1]+x_raw[3]);
          x_err[3] = (1-e1)*x_raw[3] + 0.5*e2*x_raw[2];
          
          x_err[4] = (1-e1)*x_raw[4] + 0.5*e2*x_raw[8];
          x_err[5] = (1-e1)*((1-e1)*x_raw[5] + 0.5*e2*x_raw[6]) + 0.5*e2*((1-e1)*x_raw[9] + 0.5*e2*x_raw[10]);
          x_err[6] = (1-e1)*((1-e2)*x_raw[6] + e1*(x_raw[5]+x_raw[7])) + 0.5*e2*((1-e2)*x_raw[10] + e1*(x_raw[9]+x_raw[11]));
          x_err[7] = (1-e1)*((1-e1)*x_raw[7] + 0.5*e2*x_raw[6]) + 0.5*e2*((1-e1)*x_raw[11] + 0.5*e2*x_raw[10]);
          
          x_err[8] = (1-e2)*x_raw[8] + e1*(x_raw[4] + x_raw[12]);
          x_err[9] = (1-e2)*((1-e1)*x_raw[9] + 0.5*e2*x_raw[10]) + e1*((1-e1)*x_raw[5] + 0.5*e2*x_raw[6]) + e1*((1-e1)*x_raw[13] + 0.5*e2*x_raw[14]);
          x_err[10] = (1-e2)*((1-e2)*x_raw[10] + e1*(x_raw[9]+x_raw[11])) + e1*((1-e2)*x_raw[6] + e1*(x_raw[5]+x_raw[7])) + e1*((1-e2)*x_raw[14] + e1*(x_raw[13]+x_raw[15]));
          x_err[11] = (1-e2)*((1-e1)*x_raw[11] + 0.5*e2*x_raw[10]) + e1*((1-e1)*x_raw[7] + 0.5*e2*x_raw[6]) + e1*((1-e1)*x_raw[15] + 0.5*e2*x_raw[14]);
          
          x_err[12] = (1-e1)*x_raw[12] + 0.5*e2*x_raw[8];
          x_err[13] = (1-e1)*((1-e1)*x_raw[13] + 0.5*e2*x_raw[14]) + 0.5*e2*((1-e1)*x_raw[9] + 0.5*e2*x_raw[10]);
          x_err[14] = (1-e1)*((1-e2)*x_raw[14] + e1*(x_raw[13]+x_raw[15])) + 0.5*e2*((1-e2)*x_raw[10] + e1*(x_raw[9]+x_raw[11]));
          x_err[15] = (1-e1)*((1-e1)*x_raw[15] + 0.5*e2*x_raw[14]) + 0.5*e2*((1-e1)*x_raw[11] + 0.5*e2*x_raw[10]);
          
          // save to table
          emmission_lookup[m1-1][m2-1][z][i] = x_err;
        }
        
      }
    }
  }
}

//------------------------------------------------
// MCMC::
// update transition probability lookup table. This table is a list over L-1 loci, and each element of the list is a matrix of transition probabilities of going from one HMM state to another. Note, these probabilities give the chance of moving states from the current SNP to the next SNP, which is why we have L-1 matrices for L loci (i.e. there is no matrix for the final SNP because there is nowhere to move to). This function employs the R function getTransProbs() to get eigenvalues and eigenvectors of the rate matrix, then uses these values to calculate transition probabilities for any given distance between SNPs.
void MCMC::update_transition_lookup(double f, double rho, int k, int m1, int m2, Rcpp::Function getTransProbs) {
  
  // get z_max
  int z_max = (m1 < m2) ? m1 : m2;
  
  // get eigenvalue solutions to rate matrix
  Rcpp::List Elist = getTransProbs(f, rho, k, z_max);
  vector<double> Evalues = Rcpp_to_vector_double(Elist["Evalues"]);
  vector< vector<double> > Evectors = Rcpp_to_mat_double(Elist["Evectors"]);
  vector< vector<double> > Esolve = Rcpp_to_mat_double(Elist["Esolve"]);
  
  // populate transition matrix
  for (int j=0; j<(L-1); j++) {
    // clear existing values
    for (int z1=0; z1<(m_max+1); z1++) {
      fill(transition_lookup[j][z1].begin(), transition_lookup[j][z1].end(), 0);
    }
    // populate lookup table based on rate matrix solution
    for (int z1=0; z1<(z_max+1); z1++) {
      for (int z2=0; z2<(z_max+1); z2++) {
        
        // the probability of moving from state z1 to state z2 is the sum over some weighted exponentials
        for (int i=0; i<(z_max+1); i++) {
          if (SNP_dist[j] > 0) {
            transition_lookup[j][z1][z2] += Evectors[z2][i]*Esolve[i][z1] * exp(Evalues[i] * SNP_dist[j]);
          } else {    // SNP_dist of -1 indicates jump over contigs, i.e. infinite distance
            if (Evalues[i]==0) {
              transition_lookup[j][z1][z2] += Evectors[z2][i]*Esolve[i][z1];
            }
          }
        }
        
      }
    }
  }
  
} // end of transition_update function 


//------------------------------------------------
// MCMC::
// the forward algorithm of the HMM. 
// This function does two things: 1) returns the log-likelihood integrated over all paths through the HMM, 
//                                2) updates the forward matrix. 
// The forward matrix is normalised at each step to sum to 1 (to avoid underflow issues), but this does not affect the log-likelihood calculation.
//  Remember, forward algorithm gives us the likelihood. Remember that element frwrd is defined as Pr(x_1,x_2,...,x_i,State_i=1), 

double MCMC::forward_alg(int m1, int m2) {
  
  // get z_max
  int z_max = (m1 < m2) ? m1 : m2;
  
  // clear frwrd_mat
  for (int i=0; i<int(frwrd_mat.size()); i++) {
    fill(frwrd_mat[i].begin(), frwrd_mat[i].end(), 0);
  }
  
  printf("Line 469 \n");
  // carry out first step of algorithm
  double frwrd_sum = 0;
  double logLike = 0;
  printf("Line 473 \n");
  for (int z=0; z<(z_max+1); z++) {
    frwrd_mat[z][0] = R::dbinom(z,z_max,f,false) * emmission_lookup[m1-1][m2-1][z][0][x[0]];
//    frwrd_mat[z][0] = 0.5; debug line
    frwrd_sum += frwrd_mat[z][0];
    printf("Line 469 \n");
  }
  
  logLike += log(frwrd_sum);
  for (int z=0; z<(z_max+1); z++) {
    frwrd_mat[z][0] /= frwrd_sum;
  }
  // carry out remaining steps of algorithm
  int lociiter = 1; // debug
  for (int j=1; j<L; j++) {
    frwrd_sum = 0;
    printf("Line 479 \n");
    for (int z=0; z<(z_max+1); z++) {
      // frwrd_mat[z][j] takes input from all states in iteration j-1
      for (int i=0; i<(z_max+1); i++) {
        frwrd_mat[z][j] += frwrd_mat[i][j-1] * transition_lookup[j-1][i][z];
        
      printf("Line 483 \n");
      }
      
      printf("Line 485 \n");
      printf("This is the emm prob to be multiplied");
      print((emmission_lookup[m1-1][m2-1][z][j][x[j]]));
      
      
      frwrd_mat[z][j] *= emmission_lookup[m1-1][m2-1][z][j][x[j]];
      
      
      printf("Line 487 \n");
      
      
      frwrd_sum += frwrd_mat[z][j];
      
      
      printf("Line 489 \n");
    }
    logLike += log(frwrd_sum);
    for (int z=0; z<(z_max+1); z++) {
      frwrd_mat[z][j] /= frwrd_sum;
      printf("Line 498 \n");
    }
    printf("This is iteration \n");
    print(lociiter);
    lociiter++;
  }
  
  return(logLike);
}

//------------------------------------------------
// MCMC::
// the backward algorithm of the HMM. 
// This function updates the backwards matrix, while normalising at each step to sum to 1 (to avoid underflow issues).
// Remember, The element bkwrd contains the probability of all data past observation x[i], given that we are in state 1 at time i. 
// In other words it contains Pr(x_{i+1},x_{i+1},...x_n | State_i=1).

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
      // bkwrd_mat[z][j] contains the probability of all data after this point, and therefore takes input from all states that can be reached from this one
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
// calculate IBD matrix from forward and backward matrices. This is simple product of matrices, but is normalised to sum to 1 over HMM states.
void MCMC::get_IBD() {
  
  // get z_max
  int z_max = (m1 < m2) ? m1 : m2;
  
  f_ind = 0;
  double Lcomb = 0;
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
  //    f_ind += IBD_mat[1][j]; // original when only considering MOI 1,1

    for (int z=1; z<(z_max+1); z++){
      f_ind += IBD_mat[z][j] * z * SNP_dist[j]; // AUC -- z+1 to include the zero level
    }
    
    Lcomb += z_max*SNP_dist[j]; // AUC -- z+1 to include the zero level
    
  }
 //   f_ind /= double(L);  // original when only considering MOI 1,1
  f_ind /= double(Lcomb); 

  /*
  double state_prob_sum = 0;
  vector<double> state_prob(z_max+1);
  for (int z=0; z<(z_max+1); z++) {
  state_prob[z] = R::dbinom(z,z_max,f,false) * emmission_lookup[m1-1][m2-1][z][0][x[0]];
  state_prob_sum += state_prob[z];
  }
  int sim_state = sample1(state_prob, state_prob_sum) - 1;
  
  sim_trans_n = 0;
  for (int j=1; j<L; j++) {
  state_prob_sum = 0;
  for (int z=0; z<(z_max+1); z++) {
  state_prob[z] = transition_lookup[j-1][sim_state][z];
  state_prob[z] *= emmission_lookup[m1-1][m2-1][z][j][x[j]];
  state_prob_sum += state_prob[z];
  }
  int sim_state_new = sample1(state_prob, state_prob_sum) - 1;
  if (sim_state_new != sim_state) {
  sim_trans_n++;
  }
  sim_state = sim_state_new;
  }
  */
}

//------------------------------------------------
// MCMC::
// propose new value of m given current value and sampling weights. New value cannot be 0 or greater than m_max.
double MCMC::propose_m(double m_current, double weight_move, double weight_stay) {
  
  double m_prop = m_current;
  if (rbernoulli1( weight_move/double(weight_move + weight_stay) )) {
    m_prop = rbernoulli1(0.5) ? m_current+1 : m_current-1;
    m_prop = (m_prop==0 || m_prop>m_max) ? m_current : m_prop;
  }
  
  return(m_prop);
}

//------------------------------------------------
// MCMC::
// New value cannot be 0 or greater than k_max and k is proposed given current k and sampling weights
double MCMC::propose_k(double k_current, double weight_move, double weight_stay) {
  
  double k_prop = k_current;
  if (rbernoulli1( weight_move/double(weight_move + weight_stay) )) {
    k_prop = rbernoulli1(0.5) ? k_current+1 : k_current-1;
    k_prop = (k_prop==0 || k_prop>k_max) ? k_current : k_prop;
  }
  
  return(k_prop);
}
