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
  int L = int(loci.size());
  int N = Rcpp_to_int(args["N"]);
  double rho = Rcpp_to_double(args["rho"]);
  int generations = Rcpp_to_int(args["generations"]);
  
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
  map<int, pair<int, int>> locus_map;
  for (int t=0; t<int(1e9); ++t) {
    
    // draw events in all samples and loci
    for (int i=0; i<n; ++i) {
      locus_map.clear();
      for (int l=0; l<L; ++l) {
        
        // skip if inactive
        if (!active[i][l]) {
          continue;
        }
        
        // if first time seeing this value in this sample then always draw new
        // parent. Otherwise draw new parent based on recombination probability
        // and distance since previous locus of this type
        int new_value = 0;
        if (locus_map.count(x[i][l]) != 1) {  // first time seeing this value
          
          // draw new parent
          new_value = sample2(0,N-1);
          
        } else {  // seen before
          
          // get probability of recombination
          int delta = loci[l] - locus_map[x[i][l]].first;
          double prob_break = 1 - exp(-rho*delta);
          
          // if recombination then draw new parent, otherwise same parent as
          // previous locus with this value
          if (rbernoulli1(prob_break)) {
            
            // draw new parent
            new_value = sample2(0,N-1);
            
          } else {
            
            // parent equals previous locus with this value
            new_value = locus_map[x[i][l]].second;
            
          }
          
        }
        
        // store position of this locus
        locus_map[x[i][l]] = {loci[l], new_value};
        x[i][l] = new_value;
        
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
    if (total_coalesced == (n-1)*L || t == (generations-1)) {
      break;
    }
  }
  
  // return list
  return Rcpp::List::create(Rcpp::Named("coalesced") = coalesced,
                            Rcpp::Named("coalesced_time") = coalesced_time,
                            Rcpp::Named("mutated") = mutated,
                            Rcpp::Named("mutated_time") = mutated_time);
  
}