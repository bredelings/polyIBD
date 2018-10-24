
#pragma once

//------------------------------------------------
// MCMC class
class MCMC {

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
  std::vector< std::vector<  int> > m1_mat;
  std::vector< std::vector<  int> > m2_mat;

  // MCMC parameters
  int burnin;
  int samples;
  int reportIteration;

  // lookup tables
  std::vector< std::vector< std::vector< std::vector< std::vector<double> > > > > emmission_lookup;
  std::vector< std::vector< std::vector<double> > > transition_lookup;

  // transient MCMC objects
  int k, z_max;
  double f, logLike_old;
  std::vector<int> m1_eff;
  std::vector<int> m2_eff;
  std::vector<int> z_max_eff;
  std::vector< std::vector<  double> > frwrd_mat;
  std::vector< std::vector<  double> > bkwrd_mat;
  std::vector< std::vector< double> > IBD_mat;

  // objects for storing MCMC results
  std::vector<double> logLike_burnin_store;
  std::vector<double> logLike_store;
  std::vector<double> f_store;
  std::vector<double> k_store;
  std::vector< std::vector<  double> > IBD_marginal;
  int accept_rate;

  // temp objects
  double fbs;
  std::vector<double> fbs_store;
 // int sim_trans_n;
 // std::vector<int> sim_trans_n_store;

  // misc objects
  double f_propSD;
  int k_weight_stay;
  int k_weight_move;
  int IBD_index;


  // PUBLIC FUNCTIONS

  // constructors
  MCMC(Rcpp::List args, Rcpp::List args_functions);

  // other functions
  void burnin_MCMC(Rcpp::List args_functions);
  void samp_MCMC(Rcpp::List args_functions);
  void define_emmission_lookup();
  void update_transition_lookup(double f, double rho, int k, std::vector<int> z_max_eff, Rcpp::Function getTransProbs);
  double forward_alg(std::vector<int> z_max_eff, std::vector<int> m1_eff, std::vector<int> m2_eff);
  void backward_alg(std::vector<int> z_max_eff, std::vector<int> m1_eff, std::vector<int> m2_eff);
  void get_IBD();
  double propose_k(double k_current, double weight_move, double weight_stay);
};
