#include <Rcpp.h>
#include <RcppParallel.h>
#include "SMC.h"
#include "misc.h"
#include "probability.h"

using namespace std;

//------------------------------------------------
// draw from complete ancestral recombination graph under continent-island model
// [[Rcpp::export]]
Rcpp::List continentisland_ARG_cpp(Rcpp::List args) {
  
  // get inputs
  int n = Rcpp_to_int(args["n"]);
  vector<int> loci = Rcpp_to_vector_int(args["loci"]);
  int L = int(loci.size());
  int N = Rcpp_to_int(args["N"]);
  double rho = Rcpp_to_double(args["rho"]);
  double mu = Rcpp_to_double(args["mu"]);
  int generations = Rcpp_to_int(args["generations"]);
  
  // initialise samples. The value in x[i][j] indicates the member of the
  // population that the genetic element in sample i at locus j occupies
  vector<vector<int>> x(n);
  for (int i=0; i<n; ++i) {
    x[i] = vector<int>(L, i);
  }
  
  int next_migrant = N;
  
  vector<vector<bool>> active(n, vector<bool>(L, true));
  vector<vector<bool>> coalesced(n, vector<bool>(L, false));
  vector<vector<bool>> migrated(n, vector<bool>(L, false));
  vector<vector<int>> cw(n, vector<int>(L));
  vector<vector<int>> ct(n, vector<int>(L));
  vector<vector<int>> mt(n, vector<int>(L));
  
  // a map is used to keep track of recombination events (see below)
  map<int, pair<int, int>> recom_map;
  
  // loop backwards through time
  for (int t=0; t<generations; ++t) {
    
    // loop through all samples
    for (int i=0; i<n; ++i) {
      
      // clear the recombination map for this sample
      recom_map.clear();
      
      // loop through all loci
      for (int l=0; l<L; ++l) {
        
        // skip if inactive
        if (!active[i][l]) {
          continue;
        }
        
        // if first time seeing the value x[i][l] as we scan across loci then we
        // are essentially jumping to a new host, therefore draw new parent. If
        // we have seen x[i][l] before then we have already drawn a new parent
        // for the genetic material in this host which will be stored in the
        // recom_map. Therefore, draw new or existing parent based on
        // recombination probability
        int new_value = 0;
        if (recom_map.count(x[i][l]) != 1) {  // first time seeing this value
          
          // draw new parent
          if (rbernoulli1(mu)) {
            active[i][l] = false;
            migrated[i][l] = true;
            new_value = next_migrant++;
            mt[i][l] = t;
          } else {
            new_value = sample2(0,N-1);
          }
          
        } else {  // seen before
          
          // get probability of recombination
          int delta = loci[l] - recom_map[x[i][l]].first;  // distance since previous value
          double prob_break = 1 - exp(-rho*delta);
          
          // if recombination then draw new parent, otherwise same parent as
          // previous locus with this value
          if (rbernoulli1(prob_break)) {
            if (rbernoulli1(mu)) {
              active[i][l] = false;
              migrated[i][l] = true;
              new_value = next_migrant++;
              mt[i][l] = t;
            } else {
              new_value = sample2(0,N-1);
            }
          } else {
            new_value = recom_map[x[i][l]].second;  // the stored parent
            if (migrated[i][l-1]) {
              active[i][l] = false;
              migrated[i][l] = true;
            }
          }
          
        }  // end recom draw
        
        // update recom_map and x
        recom_map[x[i][l]] = {loci[l], new_value};
        x[i][l] = new_value;
        
      }
    }
    
    // check for coalescence
    for (int i=1; i<n; ++i) {
      for (int l=0; l<L; ++l) {
        
        // skip if i inactive
        if (!active[i][l]) {
          continue;
        }
        
        // compare against all other samples
        for (int j=0; j<i; ++j) {
          
          // skip if j inactive
          if (!active[j][l]) {
            continue;
          }
          
          if (x[i][l] == x[j][l]) {
            
            // register coalescence
            active[i][l] = false;
            coalesced[i][l] = true;
            cw[i][l] = j;
            ct[i][l] = t;
          }
        }
      }
    }  // end check for coalescence
    
  } // end loop through time
  
  // return list
  return Rcpp::List::create(Rcpp::Named("coalesced") = coalesced,
                            Rcpp::Named("coalesced_with") = cw,
                            Rcpp::Named("coalesced_time") = ct,
                            Rcpp::Named("migrated") = migrated,
                            Rcpp::Named("migrated_with") = x,
                            Rcpp::Named("migrated_time") = mt);
  
}

//------------------------------------------------
// draw from conditional SMC under continent-island model
// [[Rcpp::export]]
Rcpp::List continentisland_SMC_conditional_cpp(Rcpp::List args) {
  
  // get inputs
  bool first_draw = Rcpp_to_bool(args["first_draw"]);
  vector<vector<bool>> recom = Rcpp_to_mat_bool(args["recom"]);
  vector<vector<int>> recom_time = Rcpp_to_mat_int(args["recom_time"]);
  vector<vector<int>> event_time = Rcpp_to_mat_int(args["event_time"]);
  vector<vector<bool>> coalescence = Rcpp_to_mat_bool(args["coalescence"]);
  vector<vector<int>> coalesce_with = Rcpp_to_mat_int(args["coalesce_with"]);
  vector<vector<int>> migrate_with = Rcpp_to_mat_int(args["migrate_with"]);
  vector<int> loci = Rcpp_to_vector_int(args["loci"]);
  int L = int(loci.size());
  int N = Rcpp_to_int(args["N"]);
  double rho = Rcpp_to_double(args["rho"]);
  double m = Rcpp_to_double(args["m"]);
  int next_migrant = Rcpp_to_int(args["next_migrant"]);
  
  // initialize objects for storing new events
  vector<bool> recom_new(L, false);
  vector<int> recom_time_new(L);
  vector<int> event_time_new(L);
  vector<bool> coalescence_new(L, false);
  vector<int> coalesce_with_new(L);
  vector<int> migrate_with_new(L);
  
  // get distance between loci
  vector<int> delta(L);
  for (int l=1; l<L; ++l) {
    delta[l] = loci[l] - loci[l-1];
  }
  
  // -- special case if first draw --
  
  if (first_draw) {
    
    // draw first locus from migration time distribution
    event_time_new[0] = rgeom1(m);
    
    // draw remaining loci from Markov process
    for (int l=1; l<L; ++l) {
      
      // copy values from previous locus
      event_time_new[l] = event_time_new[l-1];
      migrate_with_new[l] = migrate_with_new[l-1];
      
      // draw recombination events, leading to new draw from migration times
      double prob_recom = 1 - exp(-rho*(event_time_new[l] + 1)*delta[l]);
      recom_new[l] = rbernoulli1(prob_recom);
      if (recom_new[l]) {
        int split_time = sample2(0, event_time_new[l]);
        recom_time_new[l] = split_time;
        event_time_new[l] = split_time + rgeom1(m);
        migrate_with_new[l] = next_migrant++;
      }
    }
    
  } else {  // -- if not first draw --
    
    // misc objects
    int n = event_time[0].size();  // number of existing genotypes
    vector<int> t_sort;  // sorted event times
    
    // -- draw first locus from standard coalescent process --
    
    // draw timing of coalescence or migration
    t_sort = event_time[0];
    sort(t_sort.begin(), t_sort.end());
    int k = n;
    int new_time = rgeom1(k/double(N) + m);
    for (int i=0; i<n; ++i) {
      if (new_time > t_sort[i]) {
        k--;
        new_time = t_sort[i] + rgeom1(k/double(N) + m);
      } else {
        break;
      }
    }
    event_time_new[0] = new_time;
    double prob_coalesce = k/(k + N*m);
    coalescence_new[0] = rbernoulli1(prob_coalesce);
    
    // draw target of coalescence or migration
    if (coalescence_new[0]) {
      int tmp1 = sample2(0,k-1);
      for (int i=0; i<n; ++i) {
        if (event_time[0][i] > new_time) {
          if (tmp1 == 0) {
            coalesce_with_new[0] = i;
            break;
          }
          tmp1--;
        }
      }
    } else {
      migrate_with_new[0] = next_migrant++;
    }
    
    // --- draw remaining loci from Markov process ---
    
    // loop through remaining loci
    for (int l=1; l<L; ++l) {
      
      // copy over previous locus values
      event_time_new[l] = event_time_new[l-1];
      coalescence_new[l] = coalescence_new[l-1];
      coalesce_with_new[l] = coalesce_with_new[l-1];
      migrate_with_new[l] = migrate_with_new[l-1];
      
      // this lineage may have coalesced to a target lineage in locus l-1 that
      // underwent a recombination event between l-1 and l, making this
      // coalescence impossible. Therefore reconfigure tree as needed to account
      // for this
      if (coalescence_new[l-1]) {
        int target = coalesce_with_new[l-1];
        if (recom[l][target] && recom_time[l][target] < event_time_new[l-1]) {
          event_time_new[l] = event_time[l-1][target];
        }
      }
      
      // draw whether there is recombination event since previous locus
      double prob_recom = 1 - exp(-rho*(event_time_new[l] + 1)*delta[l]);
      recom_new[l] = rbernoulli1(prob_recom); 
      if (recom_new[l]) {
        
        // draw recombination split time uniformly over this branch
        int split_time = sample2(0, event_time_new[l]);
        recom_time_new[l] = split_time;
        
        // draw new timing of coalescence or migration
        t_sort = event_time[l];
        sort(t_sort.begin(), t_sort.end());
        int k = n;
        for (int i=0; i<n; ++i) {
          if (t_sort[i] > split_time) {
            break;
          }
          k--;
        }
        int new_time = split_time + rgeom1(k/double(N) + m);
        for (int i=0; i<n; ++i) {
          if (t_sort[i] <= split_time) {
            continue;
          }
          if (new_time > t_sort[i]) {
            k--;
            new_time = t_sort[i] + rgeom1(k/double(N) + m);
          } else {
            break;
          }
        }
        event_time_new[l] = new_time;
        double prob_coalesce = k/(k + N*m);
        coalescence_new[l] = rbernoulli1(prob_coalesce);
        
        // draw target of coalescence or migration
        if (coalescence_new[l]) {
          int tmp1 = sample2(0,k-1);
          for (int i=0; i<n; ++i) {
            if (event_time[l][i] > new_time) {
              if (tmp1 == 0) {
                coalesce_with_new[l] = i;
                break;
              }
              tmp1--;
            }
          }
        } else {
          migrate_with_new[l] = next_migrant++;
        }
        
      }  // end if recom
      
    }  // end loop over loci
    
  }  // end if/else first draw condition
  
  
  // return list
  return Rcpp::List::create(Rcpp::Named("recom_new") = recom_new,
                            Rcpp::Named("recom_time_new") = recom_time_new,
                            Rcpp::Named("event_time_new") = event_time_new,
                            Rcpp::Named("coalescence_new") = coalescence_new,
                            Rcpp::Named("coalesce_with_new") = coalesce_with_new,
                            Rcpp::Named("migrate_with_new") = migrate_with_new);
  
}
