
#pragma once

//------------------------------------------------
// MCMC class
class MCMC {
  
public:
  
  // PUBLIC OBJECTS
  
  // data and model parameters
  std::vector<int> x;
  int n;
  std::vector<std::vector<bool>> recom;
  std::vector<std::vector<int>> recom_time;
  std::vector<std::vector<int>> event_time;
  std::vector<std::vector<bool>> coalescence;
  std::vector<std::vector<int>> coalesce_target;
  std::vector<std::vector<int>> migrate_target;
  std::vector<int> loci;
  int N;
  double rho;
  double m;
  int next_migrant;
  
  // MCMC parameters
  int burnin;
  int samples;
  
  // misc parameters
  int L;
  std::vector<int> delta;
  
  
  // PUBLIC FUNCTIONS
  
  // constructors
  MCMC(Rcpp::List args);
  
  // other functions
  //void burnin_MCMC(Rcpp::List args_functions);
  
};