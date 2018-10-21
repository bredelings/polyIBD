#include <Rcpp.h>
#include <RcppParallel.h>
#include "stgIMCMC.h"
#include "misc.h"
#include "probability.h"

using namespace std;

//------------------------------------------------
// MCMC::
// constructor for MCMC class
//------------------------------------------------
// MCMC::
// constructor for MCMC class
stgIMCMC::stgIMCMC(Rcpp::List args, Rcpp::List args_functions) {

  // data and model parameters
  x = Rcpp_to_vector_int(args["x"]);
  L = x.size();
  p_vec = Rcpp_to_vector_double(args["p"]);
  rho = Rcpp_to_double(args["rho"]);
  SNP_dist = Rcpp_to_vector_int(args["SNP_dist"]);
  e1 = Rcpp_to_double(args["e1"]);
  e2 = Rcpp_to_double(args["e2"]);
  m_max = Rcpp_to_int(args["m_max"]);
  k_max = Rcpp_to_int(args["k_max"]);

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
  f = 0.01;
  k=1;
  logLike_old = 0;
  frwrd_mat = vector< vector< double> >(m_max+1, vector< double>(L));
  bkwrd_mat = vector< vector< double> >(m_max+1, vector< double>(L));
  IBD_mat = vector< vector< double> >(m_max+1, vector<double>(L));


  // objects for storing MCMC results
  logLike_burnin_store = vector< double>(burnin);
  logLike_store = vector< double>(samples);
  m1_store = vector<int>(samples);
  f_store = vector<double>(samples);
  k_store = vector<double>(samples);
  IBD_marginal = vector< vector< double> >(m_max+1, vector< double>(L));
  effMOI = vector< vector< double> >(samples, vector<double>(L));
  accept_rate = 0;

  // temp objects
  fws = 0;
  fws_store = vector<double>(samples);
  // sim_trans_n = 0;
  // sim_trans_n_store = vector<int>(samples);

  // calculate initial likelihood

  update_transition_lookup(f, rho, k, m1, args_functions["getTransProbs"]);

  logLike_old = forward_alg(m1);


  logLike_burnin_store[0] = logLike_old;

  // misc objects
  // m weights dictate the chance of proposing a new m value ("move") vs. sticking with the current value ("stay"). These weights are updated during the burn-in phase, converging to an efficient proposal distribution. The chance of staying is never allowed to exceed 100* the chance of moving, to ensure there is always some chance of exploring new m values.
  // f_propSD and k_propSD are updated during the burn-in phase using Robbins-Monro (although k_propSD isn't being used right now).
  m1_weight_stay = 1;
  m1_weight_move = 1;
  f_propSD = 0.2;
  k_weight_stay = 1;
  k_weight_move = 1;
}



//------------------------------------------------
// stgIMCMC::
// Run burn-in phase of MCMC. Each iteration we alternate which parameters are updated.

void stgIMCMC::burnin_MCMC(Rcpp::List args_functions) {
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

    // propose m1 (MOI)
    int m1_prop = propose_m(m1, m1_weight_move, m1_weight_stay);

    // if no change in m then propose either f or k
    double f_prop = f;
	  int k_prop = k;
    if (m1_prop==m1) {
      if (rbernoulli1(0.5)) {
        f_prop = rnorm1_interval(f, f_propSD, 0, 1);
      } else {
         k_prop = propose_k(k, k_weight_move, k_weight_stay);
      }
    }

    // update transition probabilities and calculate new likelihood
    update_transition_lookup(f_prop, rho, k_prop, m1_prop, args_functions["getTransProbs"]);
    double logLike_new = forward_alg(m1_prop);


    // Metropolis-Hastings step
    // note that all proposal distributions are symmetric, therefore no Hastings step is required -- we meet g(x | x') = g(x' | x)
    // if accept
    if (log(runif_0_1()) < (logLike_new-logLike_old)) {

      // update m1 sampling weights
      m1_weight_move = (m1==m1_prop) ? m1_weight_move : ++m1_weight_move;

      // or update f_propSD
      if (m1==m1_prop && f_prop!=f) {
        f_propSD  += (1-0.23)/sqrt(double(rep));
      }
      // or update k_propSD
      if (m1==m1_prop && k_prop!=k) {
        k_weight_move = (k==k_prop) ? k_weight_move : ++k_weight_move;
      }

      // update parameter values and likelihood
      f = f_prop;
      k = k_prop;
      m1 = m1_prop;
      logLike_old = logLike_new;

    }
    // if reject
    else {

      // update m1 sampling weights
      m1_weight_stay = (m1==m1_prop) ? m1_weight_stay : ++m1_weight_stay;

      // limit m1 sampling weights
      m1_weight_stay = (m1_weight_stay > 100*m1_weight_move) ? 100*m1_weight_move : m1_weight_stay;

      // update f_propSD
      if (m1==m1_prop && f_prop!=f) {
        f_propSD  -= 0.23/sqrt(double(rep));
        f_propSD = (f_propSD < 0) ? -f_propSD : f_propSD;
      }

      // update k move
      if (m1==m1_prop && k_prop!=k) {
        k_weight_stay = (k==k_prop) ? k_weight_stay : ++k_weight_stay;
        // limit k sampling weights
        k_weight_stay = (k_weight_stay > 100*k_weight_move) ? 100*k_weight_move : k_weight_stay;
      }
    }
    // store logLike
    logLike_burnin_store[rep] = logLike_old;

  }   // end MCMC burninloop

}


//------------------------------------------------
// MCMC::
// Run sampling phase of MCMC. Most of this function is very similar to the burn-in phase.
// Differences are that sampling weights and standard deviations are not updated by Robbins-Monro,
// and a greater number of outputs are stored.
void stgIMCMC::samp_MCMC(Rcpp::List args_functions) {

  // announce burn-in phase
  print("   sampling phase");

  // calculate initial IBD matrix
  backward_alg(m1);
  get_IBD();

  // loop through sampling iterations
  for (int rep=0; rep<samples; rep++) {

    // report progress
    if (reportIteration>0) {
      if (((rep+1) % reportIteration)==0) {
        print("      iteration",rep+1);
      }
    }

    // propose m1 (MOI)
    int m1_prop = propose_m(m1, m1_weight_move, m1_weight_stay);

    // if no change in m, therefore propose either f or k
    double f_prop = f;
    int k_prop = k;
    if (m1_prop==m1) {
      if (rbernoulli1(0.5)) {
        f_prop = rnorm1_interval(f, f_propSD, 0, 1);
      } else {
        k_prop = propose_k(k, k_weight_move, k_weight_stay);
      }
    }

    // update transition probabilities and calculate new likelihood
    update_transition_lookup(f_prop, rho, k_prop, m1_prop, args_functions["getTransProbs"]);
    double logLike_new = forward_alg(m1_prop);

    // Metropolis-Hastings step
    // note that all proposal distributions are symmetric, therefore no Hastings step is required
    // if accept
    if (log(runif_0_1()) < (logLike_new-logLike_old)) {

      // update parameter values and likelihood
      f = f_prop;
      k = k_prop;
      m1 = m1_prop;
      logLike_old = logLike_new;

      // update IBD_mat
      backward_alg(m1);
      get_IBD();

      // update acceptance rate
      accept_rate ++;
    }

    // store current values
    logLike_store[rep] = logLike_old;
    m1_store[rep] = m1;
    f_store[rep] = f;
    fws_store[rep] = fws;
    k_store[rep] = k;
    //sim_trans_n_store[rep] = sim_trans_n;

    // add current IBD_mat to running estimate
    for (int i=0; i<(m_max+1); i++) {
      for (int j=0; j<L; j++) {
        IBD_marginal[i][j] += IBD_mat[i][j];
      }
    }

    // determine the effective MOI for this report
    // defined z_max in forward, backward, and get_IBD as m1-1
    // m1 - zlvl to get the effective MOI estimate
    z_maxvec = vector< double>(m1);

    for(int j=0; j<L; j++){
      for(int i=0; i<m1; i++){
        fill(z_maxvec.begin(), z_maxvec.end(), 0);
        z_maxvec[i] = IBD_mat[i][j];
      }
      effMOI[rep][j] = m1 - sampleZ(z_maxvec);
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
// stgIMCMC::
// define emmission probability lookup table.
 // This table is essentially a list over m1, then h, where h
// is the number of HMM states at each locus (i.e. h = 1 + the
// maximum level of IBD between samples). At the final level of
// the list is a matrix with L rows and 4 columns giving emmission probabilities for each locus.
// The 4 columns correspond to the 4 possible genotypes combinations (NA means missing data):
// 0 = {NA}
// 1 = {A}
// 2 = {Aa}
// 3 = {a}

void stgIMCMC::define_emmission_lookup() {

  // temporary vector for storing emmission probs
  vector<double> x_raw(4);
  vector<double> x_err(4);

  // loop through m1
  emmission_lookup = vector< vector< vector< vector<double> > > > (m_max);
  for (int m1=1; m1<=m_max; m1++) {
      // define h as 1 + maximum IBD
      int h = m1;

      // loop through z
      emmission_lookup[m1-1] = vector< vector< vector<double> > >(h);
      for (int z=0; z<h; z++) {
        emmission_lookup[m1-1][z] = vector< vector<double> >(L, vector<double>(4));

        // loop through loci
        for (int i=0; i<L; i++) {
          double p = p_vec[i];
          double q = 1-p;

          // calculate raw emmission probability (without error)
          // if no IBD
          if (z==0) {

            x_raw[0] = 1;
            x_raw[1] = pow(p,m1);
            x_raw[2] = (1 - pow(p,m1) - pow(q,m1));
            x_raw[3] = pow(q,m1);

          }
          // if some IBD
          else {

            x_raw[0] = 1;
            x_raw[1] = pow(p,m1-z);
            x_raw[2] = (1 - pow(p,m1-z) - pow(q,m1-z));
            x_raw[3] = pow(q,m1-z);
          }

          // incorporate error
          // e1 = prob homo incorrectly called as het
          // e2 = prob het incorrectly called as homo

          x_err[0] = x_raw[0];
          x_err[1] = (1-e1)*x_raw[1] + 0.5*e2*x_raw[2];
          x_err[2] = (1-e2)*x_raw[2] + e1*(x_raw[1]+x_raw[3]);
          x_err[3] = (1-e1)*x_raw[3] + 0.5*e2*x_raw[2];

          // save to table
          emmission_lookup[m1-1][z][i] = x_err;

      }
    }
  }
}


//------------------------------------------------
// stgIMCMC::
// update transition probability lookup table. This table is a list over L-1 loci, and each
// element of the list is a matrix of transition probabilities of going from one HMM state to another.
// Note, these probabilities give the chance of moving states from the current SNP to the next SNP,
// which is why we have L-1 matrices for L loci (i.e. there is no matrix for the final SNP because there
// is nowhere to move to). This function employs the R function getTransProbs() to get eigenvalues and
// eigenvectors of the rate matrix, then uses these values to calculate transition probabilities for any
// given distance between SNPs.

void stgIMCMC::update_transition_lookup(double f, double rho, int k, int m1, Rcpp::Function getTransProbs) {


  // get z_max within
  // remember because we are within the sample, there is really m1-1
  // potential pairings (assuming poisson ind processes)

  int z_max = m1-1;


  // however, this now poses the problem where if we have an MOI of 1 then there is no comparison
  // of I and U to be made. Really, a sample is fully I with itself -- this is consistent with the
  // previous notation of fws where fws > 0.95 has been used to indicate a monoclonal infxn


  if(z_max == 0){
    // populate transition matrix for this condition
    for (int j=0; j<(L-1); j++) {

      // clear existing values
      for (int z1=0; z1<(m_max); z1++) {
        fill(transition_lookup[j][z1].begin(), transition_lookup[j][z1].end(), 1); // first fill with 0s
        transition_lookup[j][1][1] = 1; // a11 is 1 becuase a strain is identical to itself
      }
    }

  } else {



    // get eigenvalue solutions to rate matrix
    Rcpp::List Elist = getTransProbs(f, rho, k, z_max);
    vector<double> Evalues = Rcpp_to_vector_double(Elist["Evalues"]);
    vector< vector<double> > Evectors = Rcpp_to_mat_double(Elist["Evectors"]);
    vector< vector<double> > Esolve = Rcpp_to_mat_double(Elist["Esolve"]);


    // populate transition matrix
    for (int j=0; j<(L-1); j++) {

      // clear existing values
      for (int z1=0; z1<(m_max); z1++) {
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
} // end of if-else conditional

//------------------------------------------------
// stgIMCMC::
// the forward algorithm of the HMM.
// This function does two things: 1) returns the log-likelihood integrated over all paths through the HMM,
//                                2) updates the forward matrix.
// The forward matrix is normalised at each step to sum to 1 (to avoid underflow issues), but this does not affect the log-likelihood calculation.
//  Remember, forward algorithm gives us the likelihood. Remember that element frwrd is defined as Pr(x_1,x_2,...,x_i,State_i=1),

double stgIMCMC::forward_alg(int m1) {

  // get z_max within
  // remember because we are within the sample, there is really m1-1
  // potential pairings (assuming poisson ind processes)

  int z_max = m1-1;

  // clear frwrd_mat
  for (int i=0; i<int(frwrd_mat.size()); i++) {
    fill(frwrd_mat[i].begin(), frwrd_mat[i].end(), 0);
  }

  // carry out first step of algorithm
  double frwrd_sum = 0;
  double logLike = 0;
  for (int z=0; z<(z_max+1); z++) {
    frwrd_mat[z][0] = R::dbinom(z,z_max,f,false) * emmission_lookup[m1-1][z][0][x[0]];
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
      // frwrd_mat[z][j] takes input from all states in iteration j-1
      for (int i=0; i<(z_max+1); i++) {
        frwrd_mat[z][j] += frwrd_mat[i][j-1] * transition_lookup[j-1][i][z];
      }

      frwrd_mat[z][j] *= emmission_lookup[m1-1][z][j][x[j]];
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
// stgIMCMC::
// the backward algorithm of the HMM.
// This function updates the backwards matrix, while normalising at each step to sum to 1 (to avoid underflow issues).
// Remember, The element bkwrd contains the probability of all data past observation x[i], given that we are in state 1 at time i.
// In other words it contains Pr(x_{i+1},x_{i+1},...x_n | State_i=1).

void stgIMCMC::backward_alg(int m1) {

  // get z_max within
  // remember because we are within the sample, there is really m1-1
  // potential pairings (assuming poisson ind processes)

  int z_max = m1-1;

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
        bkwrd_mat[z][j] += transition_lookup[j][z][i]*bkwrd_mat[i][j+1]*emmission_lookup[m1-1][i][j+1][x[j+1]];
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
void stgIMCMC::get_IBD() {

  // get z_max within
  // remember because we are within the sample, there is really m1-1
  // potential pairings (assuming poisson ind processes)

  int z_max = m1-1;

    // initialise values
    fws = 0;
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

      for (int z=1; z<(z_max+1); z++){
        fws += IBD_mat[z][j] * z * SNP_dist[j]; // AUC -- z+1 to include the zero level
      }

      Lcomb += z_max*SNP_dist[j]; // AUC -- z+1 to include the zero level

    }

    if(Lcomb != 0){
      fws /= double(Lcomb);
    }
  }


//------------------------------------------------
// stgIMCMC::
// propose new value of m given current value and sampling weights. New value cannot be 0 or greater than m_max.
double stgIMCMC::propose_m(double m_current, double weight_move, double weight_stay) {
  double m_prop = m_current;
  if (rbernoulli1( weight_move/double(weight_move + weight_stay) )) {
    m_prop = rbernoulli1(0.5) ? m_current+1 : m_current-1;
    m_prop = (m_prop==0 || m_prop>m_max) ? m_current : m_prop;
  }
  return(m_prop);
}


//------------------------------------------------
// stgIMCMC::
// propose new value of K given current value and sampling weights. New value cannot be 0 or greater than m_max.
double stgIMCMC::propose_k(double k_current, double weight_move, double weight_stay) {
  double k_prop = k_current;
  if (rbernoulli1( weight_move/double(weight_move + weight_stay) )) {
    k_prop = rbernoulli1(0.5) ? k_current+1 : k_current-1;
    k_prop = (k_prop==0 || k_prop>k_max) ? k_current : k_prop;
  }
  return(k_prop);
}
