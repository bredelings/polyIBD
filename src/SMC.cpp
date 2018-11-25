
#include <Rcpp.h>
#include <RcppParallel.h>
#include "SMC.h"
#include "misc.h"
#include "probability.h"

using namespace std;

//------------------------------------------------
// draw from complete ancestral recombination graph under wright-fisher model
// [[Rcpp::export]]
Rcpp::List wrightfisher_ARG_cpp(Rcpp::List args) {
  
  // get inputs
  int n = Rcpp_to_int(args["n"]);
  vector<int> loci = Rcpp_to_vector_int(args["loci"]);
  int N = Rcpp_to_int(args["N"]);
  double rho = Rcpp_to_double(args["rho"]);
  int generations = Rcpp_to_int(args["generations"]);
  bool generations_infinite = (generations == 0);
  
  // get genetic distances
  int L = int(loci.size());
  vector<int> delta(L);
  for (int l=1; l<L; ++l) {
    delta[l] = loci[l]-loci[l-1];
  }
  
  // initialise samples. The value in x[i][j] indicates the member of the
  // population that sample i at locus j occupies
  vector<vector<int>> x(n);
  for (int i=0; i<n; ++i) {
    x[i] = vector<int>(L, i);
  }
  
  // store whether sample has coalesced or mutated at each locus, plus timing of
  // these events
  vector<vector<bool>> coalesced(n, vector<bool>(L, false));
  vector<vector<bool>> mutated(n, vector<bool>(L, false));
  vector<vector<int>> coalesced_time(n, vector<int>(L));
  vector<vector<int>> mutated_time(n, vector<int>(L));
  int total_coalesced = 0;
  
  // mask indicating whether genetic element is still active (i.e. un-coalesced
  // and un-mutated)
  vector<vector<bool>> active(n, vector<bool>(L, true));
  
  // loop backwards through time
  int t = -1;
  while (true) {
    t++;
    
    // draw events in all samples and loci
    for (int i=0; i<n; ++i) {
      for (int l=0; l<L; ++l) {
        
        // skip if inactive
        if (!active[i][l]) {
          continue;
        }
        
        // draw whether recombination breakpoint since previous locus
        double prob_break = 1 - exp(-rho*delta[l]);
        if (l == 0) {
          prob_break = 1;
        }
        if (rbernoulli1(prob_break)) {
          
          // draw new parent
          x[i][l] = sample2(0,N-1);
          
        } else {
          
          // parent equals previous locus
          x[i][l] = x[i][l-1];
          
        }
      }
    }
    
    // check for coalescence
    for (int i=1; i<n; ++i) {
      for (int l=0; l<L; ++l) {
        
        // skip if inactive
        if (!active[i][l]) {
          continue;
        }
        
        // compare against all other samples
        for (int j=0; j<i; ++j) {
          if (x[i][l] == x[j][l]) {
            
            // register coalescence
            coalesced[i][l] = true;
            coalesced_time[i][l] = t;
            active[i][l] = false;
            total_coalesced++;
            
          }
        }
      }
    }
    
    // break if all coalesced or if reached max time
    if (total_coalesced == (n-1)*L || (!generations_infinite && t == (generations-1))) {
      break;
    }
  }
  
  // return list
  return Rcpp::List::create(Rcpp::Named("coalesced") = coalesced,
                            Rcpp::Named("coalesced_time") = coalesced_time,
                            Rcpp::Named("mutated") = mutated,
                            Rcpp::Named("mutated_time") = mutated_time);
  
}