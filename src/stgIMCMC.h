
#pragma once

//------------------------------------------------
// MCMC class
class stgIMCMC {

public:

  // PUBLIC OBJECTS

  // data and model parameters
  std::vector<int> x;
  int L;
  double rho;
  std::vector<double> p_vec;
  std::vector<int> SNP_dist;
  double e1;
  double e2;
  int m_max;
  int k_max;

  // MCMC parameters
  int burnin;
  int samples;
  int reportIteration;

  // lookup tables
  std::vector< std::vector< std::vector< std::vector<double> > > > emmission_lookup;
  std::vector< std::vector< std::vector<double> > > transition_lookup;

  // transient MCMC objects
  int m1, m2, k, z_max;
  double f, logLike_old;
  double fws;
  // int sim_trans_n;
  std::vector< std::vector<  double> > frwrd_mat;
  std::vector< std::vector<  double> > bkwrd_mat;
  std::vector< std::vector< double> > IBD_mat;

  // objects for storing MCMC results
  std::vector< double> logLike_burnin_store;
  std::vector< double> logLike_store;
  std::vector<int> m1_store;
  std::vector<double> f_store;
  std::vector<double> k_store;
  std::vector< std::vector<  double> > IBD_marginal;
  int accept_rate;
  std::vector<double> fws_store;
 // std::vector<int> sim_trans_n_store;

  // misc objects
  int m1_weight_stay;
  int m1_weight_move;
  double f_propSD;
  int k_propSD;
  int IBD_index;


  // PUBLIC FUNCTIONS

  // constructors
  stgIMCMC(Rcpp::List args, Rcpp::List args_functions);

  // other functions
  void burnin_MCMC(Rcpp::List args_functions);
  void samp_MCMC(Rcpp::List args_functions);
  void define_emmission_lookup();
  void update_transition_lookup(double f, double rho, int k, int m1, Rcpp::Function getTransProbs);
  double forward_alg(int m1);
  void backward_alg(int m1);
  void get_IBD();
  double propose_m(double m_current, double weight_move, double weight_stay);
  double propose_k(double k_current, double weight_move, double weight_stay);
};
