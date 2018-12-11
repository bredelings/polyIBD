
#include <Rcpp.h>
#include "MCMC.h"
#include "misc.h"
#include "probability.h"

using namespace std;

//------------------------------------------------
// constructor
MCMC::MCMC(Rcpp::List args) {
  
  // data and model parameters
  x = rcpp_to_vector_int(args["x"]);
  n = rcpp_to_int(args["n"]);
  recom = rcpp_to_matrix_bool(args["recom"]);
  recom_time = rcpp_to_matrix_int(args["recom_time"]);
  event_time = rcpp_to_matrix_int(args["event_time"]);
  coalescence = rcpp_to_matrix_bool(args["coalescence"]);
  coalesce_target = rcpp_to_matrix_int(args["coalesce_target"]);
  migrate_target = rcpp_to_matrix_int(args["migrate_target"]);
  loci = rcpp_to_vector_int(args["loci"]);
  N = rcpp_to_int(args["N"]);
  rho = rcpp_to_double(args["rho"]);
  m = rcpp_to_double(args["m"]);
  next_migrant = rcpp_to_int(args["next_migrant"]);
  
  // MCMC parameters
  burnin = rcpp_to_int(args["burnin"]);;
  samples = rcpp_to_int(args["samples"]);;
  
  // misc parameters
  L = int(loci.size());
  delta = vector<int>(L);
  
  // get distance between loci
  for (int l=1; l<L; ++l) {
    delta[l] = loci[l] - loci[l-1];
  }
  
}