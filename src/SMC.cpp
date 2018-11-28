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
  vector<vector<int>> coalesced_time(n, vector<int>(L));
  int total_coalesced = 0;
  
  // a map is used to keep track of unique values in each sample
  map<int, pair<int, int>> locus_map;
  
  // loop backwards through time
  for (int t=0; t<int(1e9); ++t) {
    
    // draw events in all samples and loci
    for (int i=0; i<n; ++i) {
      locus_map.clear();
      for (int l=0; l<L; ++l) {
        
        // skip if already coalesced
        if (coalesced[i][l]) {
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
            new_value = sample2(0,N-1);
          } else {
            new_value = locus_map[x[i][l]].second;
          }
          
        }
        
        // store new value and latest position
        locus_map[x[i][l]] = {loci[l], new_value};
        x[i][l] = new_value;
        
      }
    }
    
    // check for coalescence
    for (int i=1; i<n; ++i) {
      for (int l=0; l<L; ++l) {
        
        // skip if already coalesced
        if (coalesced[i][l]) {
          continue;
        }
        
        // compare against all other samples
        for (int j=0; j<i; ++j) {
          if (x[i][l] == x[j][l]) {
            
            // register coalescence
            coalesced[i][l] = true;
            coalesced_time[i][l] = t;
            if (j == 0) {
              coalesced[j][l] = true;
              coalesced_time[j][l] = t;
            }
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
                            Rcpp::Named("coalesced_time") = coalesced_time);
  
}

//------------------------------------------------
// draw from sequentially Markov coalescent under wright-fisher model
// [[Rcpp::export]]
Rcpp::List wrightfisher_SMC_cpp(Rcpp::List args) {
  
  // get inputs
  int n = Rcpp_to_int(args["n"]);
  vector<int> loci = Rcpp_to_vector_int(args["loci"]);
  int L = int(loci.size());
  int N = Rcpp_to_int(args["N"]);
  double rho = Rcpp_to_double(args["rho"]);
  
  // initialise samples. The value in x[l][i] stores the time at which sample i
  // coalesced with one of the other samples at locus l.
  vector<vector<int>> cw(L, vector<int>(n, -1));
  vector<vector<int>> ct(L, vector<int>(n));
  
  // the value T[l] stores the total branch length at locus l
  vector<int> T(L);
  
  // draw first locus from standard coalescent
  int t = 0;
  for (int i=0; i<(n-1); ++i) {
    int k = (n-i);
    t += rgeom1(0.5*k*(k-1)/double(N));
    cw[0][i] = sample2(i+1,n-1);
    ct[0][i] = t;
  }
  ct[0][n-1] = ct[0][n-2];
  T[0] = sum(ct[0]);
  
  // loop through remaining loci
  for (int l=1; l<L; ++l) {
    
    // initially make tree same as previous locus
    cw[l] = cw[l-1];
    ct[l] = ct[l-1];
    T[l] = T[l-1];
    
    // draw recombination events until reach next locus
    double p = loci[l-1] + rexp1(rho*T[l]);
    while (p < loci[l]) {
      
      // draw timing of recombination event uniformally over branch lengths
      int recom_lineage = sample1(ct[l], T[l])-1;
      int recom_t1 = sample2(0,ct[l][recom_lineage]);
      
      // find the first lineage to coalesce with recom_lineage above the point
      // recom_t1. There may be no such lineage, in which case we can make
      // recom_lineage "float" without any further modifications.
      int first_lineage = -1;
      int first_lineage_time = 0;
      for (int i=0; i<n; ++i) {
        if (i == recom_lineage) {
          continue;
        }
        if (cw[l][i] == recom_lineage && ct[l][i] >= recom_t1) {
          if (ct[l][i] < first_lineage_time || first_lineage == -1) {
            first_lineage = i;
            first_lineage_time = ct[l][i];
          }
        }
      }
      
      // If there is a first_lineage then it must "take over" from recom_lineage
      // in all coalescences above recom_t1.
      if (first_lineage != -1) {
        cw[l][first_lineage] = cw[l][recom_lineage];
        ct[l][first_lineage] = ct[l][recom_lineage];
        for (int i=0; i<n; ++i) {
          if (cw[l][i] == recom_lineage) {
            cw[l][i] = first_lineage;
          }
        }
      }
      
      // mask out recom_lineage
      cw[l][recom_lineage] = -1;
      ct[l][recom_lineage] = -1;
      
      // find number of extant lineages at time recom_t1
      vector<int> ct_increasing = ct[l];
      sort(ct_increasing.begin(), ct_increasing.end());
      int k = n;
      for (int i = 1; i<n; ++i) {
        k--;
        if (ct_increasing[i] > recom_t1) {
          break;
        }
      }
      
      // draw time at which floating lineage reconnects with tree
      int recom_t2 = recom_t1 + rgeom1(k/double(N));
      while (recom_t2 >= ct_increasing[n-k] && k > 1) {
        recom_t2 = ct_increasing[n-k] + rgeom1((k-1)/double(N));
        k--;
      }
      
      // draw lineage that floating lineage connects to
      int float_connect = 0;
      int rand1 = sample2(0,k-1);
      for (int i = 0; i<n; ++i) {
        if (i == recom_lineage) {
          continue;
        }
        if (ct[l][i] >= recom_t2 || cw[l][i] == -1) {
          if (rand1 == 0) {
            float_connect = i;
          }
          rand1--;
        }
      }
      
      // reconnect recom_lineage with tree
      cw[l][recom_lineage] = float_connect;
      ct[l][recom_lineage] = recom_t2;
      
      // update the final coalescent time as needed
      int trunk = 0;
      int coal_final = 0;
      for (int i = 0; i<n; ++i) {
        if (cw[l][i] == -1) {
          trunk = i;
        } else {
          if (ct[l][i] > coal_final) {
            coal_final = ct[l][i];
          }
        }
      }
      ct[l][trunk] = coal_final;
      T[l] = sum(ct[l]);
      
      // draw next recombination position
      p += rexp1(rho*T[l]);
      
    }  // end while loop
    
  }  // end loop through loci
  
  // return list
  return Rcpp::List::create(Rcpp::Named("coalesced_time") = ct);
  
}