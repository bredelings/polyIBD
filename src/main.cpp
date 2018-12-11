
#include <chrono>
#include "main.h"
#include "misc_v2.h"
#include "probability.h"
#include "MCMC.h"

using namespace std;

//------------------------------------------------
// draw from complete ancestral recombination graph under continent-island model
// [[Rcpp::export]]
Rcpp::List continentisland_ARG_cpp(Rcpp::List args) {
  
  // get inputs
  int n = rcpp_to_int(args["n"]);
  vector<int> loci = rcpp_to_vector_int(args["loci"]);
  int L = int(loci.size());
  int N = rcpp_to_int(args["N"]);
  double rho = rcpp_to_double(args["rho"]);
  double m = rcpp_to_double(args["m"]);
  int generations = rcpp_to_int(args["generations"]);
  
  // initialise samples. The value in x[i][j] indicates the member of the
  // population that the genetic element in sample i at locus j occupies
  vector<vector<int>> x(n);
  for (int i=0; i<n; ++i) {
    x[i] = vector<int>(L, i);
  }
  
  // give migrants unique indices, starting with one more than the possible
  // index of individuals within the subpopulation
  int next_migrant = N;
  
  // initialise matrices for storing events and timings
  vector<vector<bool>> recom(n, vector<bool>(L, false));
  vector<vector<int>> recom_time(n, vector<int>(L));
  vector<vector<bool>> active(n, vector<bool>(L, true));
  vector<vector<bool>> coalescence(n, vector<bool>(L, false));
  vector<vector<int>> coalesce_target(n, vector<int>(L));
  vector<vector<int>> migrate_target(n, vector<int>(L));
  vector<vector<int>> event_time(n, vector<int>(L));
  
  // a map is used to keep track of recombination events (see implementation below)
  map<int, pair<int, int>> recom_map;
  
  // loop backwards through time
  for (int t=0; t<generations; ++t) {
    
    // loop through all samples
    for (int i=0; i<n; ++i) {
      
      // clear the recombination map
      recom_map.clear();
      
      // loop through all loci
      for (int l=0; l<L; ++l) {
        
        // skip if inactive
        if (!active[i][l]) {
          continue;
        }
        
        // if first time seeing the value x[i][l] then we are essentially
        // jumping to a new host, therefore draw new parent. If we have seen
        // x[i][l] before then we have already drawn a new parent for this host
        // which will be stored in the recom_map, in which case either draw new
        // parent or copy existing parent based on recombination probability
        if (recom_map.count(x[i][l]) != 1) {  // if first time seeing this value
          
          // draw new parent from within the subpopulation or from the migrant
          // pool
          if (rbernoulli1(m)) {
            event_time[i][l] = t;
            active[i][l] = false;
            migrate_target[i][l] = next_migrant;
            
            recom_map[x[i][l]] = {loci[l], next_migrant};
            x[i][l] = next_migrant++;
          } else {
            int new_parent = sample2(0,N-1);
            recom_map[x[i][l]] = {loci[l], new_parent};
            x[i][l] = new_parent;
          }
          
        } else {  // if seen this value before
          
          // get probability of recombination
          int delta = loci[l] - recom_map[x[i][l]].first;  // distance since previous value
          double prob_break = 1 - exp(-rho*delta);
          
          // draw whether there is recombination event
          int recom_true = rbernoulli1(prob_break);
          if (recom_true && !recom[i][l]) {
            recom[i][l] = true;
            recom_time[i][l] = t;
          }
          
          // if recombination then draw new parent, otherwise same parent as
          // stored in recombination map
          if (recom_true) {
            
            // draw new parent from within the subpopulation or from the migrant
            // pool
            if (rbernoulli1(m)) {
              event_time[i][l] = t;
              active[i][l] = false;
              migrate_target[i][l] = next_migrant;
              
              recom_map[x[i][l]] = {loci[l], next_migrant};
              x[i][l] = next_migrant++;
            } else {
              int new_parent = sample2(0,N-1);
              recom_map[x[i][l]] = {loci[l], new_parent};
              x[i][l] = new_parent;
            }
            
          } else {
            
            // get parent from recombination map
            int new_parent = recom_map[x[i][l]].second;  // the stored parent
            if (!active[i][l-1] && !coalescence[i][l-1]) {
              active[i][l] = false;
              event_time[i][l] = event_time[i][l-1];
              migrate_target[i][l] = migrate_target[i][l-1];
            }
            recom_map[x[i][l]] = {loci[l], new_parent};
            x[i][l] = new_parent;
          }
          
        }  // end recom draw
      }  // end loop through loci
    }  // end loop through samples
    
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
          
          // register coalescence
          if (x[i][l] == x[j][l]) {
            active[i][l] = false;
            coalescence[i][l] = true;
            coalesce_target[i][l] = j;
            event_time[i][l] = t;
          }
        }
      }
    }  // end check for coalescence
    
  } // end loop through time
  
  // return list
  return Rcpp::List::create(Rcpp::Named("recom") = recom,
                            Rcpp::Named("recom_time") = recom_time,
                            Rcpp::Named("event_time") = event_time,
                            Rcpp::Named("active") = active,
                            Rcpp::Named("coalescence") = coalescence,
                            Rcpp::Named("coalesce_target") = coalesce_target,
                            Rcpp::Named("migrate_target") = migrate_target);
  
}

//------------------------------------------------
// draw from conditional SMC under continent-island model
// [[Rcpp::export]]
Rcpp::List continentisland_SMC_conditional_cpp(Rcpp::List args) {
  
  // get inputs
  bool first_draw = rcpp_to_bool(args["first_draw"]);
  vector<vector<bool>> recom = rcpp_to_matrix_bool(args["recom"]);
  vector<vector<int>> recom_time = rcpp_to_matrix_int(args["recom_time"]);
  vector<vector<int>> event_time = rcpp_to_matrix_int(args["event_time"]);
  vector<vector<bool>> coalescence = rcpp_to_matrix_bool(args["coalescence"]);
  vector<vector<int>> coalesce_target = rcpp_to_matrix_int(args["coalesce_target"]);
  vector<vector<int>> migrate_target = rcpp_to_matrix_int(args["migrate_target"]);
  vector<int> loci = rcpp_to_vector_int(args["loci"]);
  int L = int(loci.size());
  int N = rcpp_to_int(args["N"]);
  double rho = rcpp_to_double(args["rho"]);
  double m = rcpp_to_double(args["m"]);
  int next_migrant = rcpp_to_int(args["next_migrant"]);
  
  // initialize objects for storing new events
  vector<bool> recom_new(L, false);
  vector<int> recom_time_new(L);
  vector<int> event_time_new(L);
  vector<bool> coalescence_new(L, false);
  vector<int> coalesce_target_new(L);
  vector<int> migrate_target_new(L);
  
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
      migrate_target_new[l] = migrate_target_new[l-1];
      
      // draw recombination events, leading to new draw from migration times
      double prob_recom = 1 - exp(-rho*(event_time_new[l] + 1)*delta[l]);
      recom_new[l] = rbernoulli1(prob_recom);
      if (recom_new[l]) {
        int split_time = sample2(0, event_time_new[l]);
        recom_time_new[l] = split_time;
        event_time_new[l] = split_time + rgeom1(m);
        migrate_target_new[l] = next_migrant++;
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
            coalesce_target_new[0] = i;
            break;
          }
          tmp1--;
        }
      }
    } else {
      migrate_target_new[0] = next_migrant++;
    }
    
    // --- draw remaining loci from Markov process ---
    
    // loop through remaining loci
    for (int l=1; l<L; ++l) {
      
      // copy over previous locus values
      event_time_new[l] = event_time_new[l-1];
      coalescence_new[l] = coalescence_new[l-1];
      coalesce_target_new[l] = coalesce_target_new[l-1];
      migrate_target_new[l] = migrate_target_new[l-1];
      
      // this lineage may have coalesced to a target lineage in locus l-1 that
      // underwent a recombination event between l-1 and l, making this
      // coalescence impossible. Therefore reconfigure tree as needed to account
      // for this
      if (coalescence_new[l-1]) {
        int target = coalesce_target_new[l-1];
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
                coalesce_target_new[l] = i;
                break;
              }
              tmp1--;
            }
          }
        } else {
          migrate_target_new[l] = next_migrant++;
        }
        
      }  // end if recom
      
    }  // end loop over loci
    
  }  // end if/else first draw condition
  
  // return list
  return Rcpp::List::create(Rcpp::Named("recom_new") = recom_new,
                            Rcpp::Named("recom_time_new") = recom_time_new,
                            Rcpp::Named("event_time_new") = event_time_new,
                            Rcpp::Named("coalescence_new") = coalescence_new,
                            Rcpp::Named("coalesce_target_new") = coalesce_target_new,
                            Rcpp::Named("migrate_target_new") = migrate_target_new);
  
}

//------------------------------------------------
// MCMC under continent-island model
// [[Rcpp::export]]
Rcpp::List SMC_MCMC_cpp(Rcpp::List args) {
  
  // create MCMC object
  MCMC mcmc(args);
  
  
  /*
  
  // initialize objects for storing new events
  vector<bool> recom_new(L, false);
  vector<int> recom_time_new(L);
  vector<int> event_time_new(L);
  vector<bool> coalescence_new(L, false);
  vector<int> coalesce_target_new(L);
  vector<int> migrate_target_new(L);
  
  
  
  // -- special case if first draw --
  
  if (first_draw) {
    
    // draw first locus from migration time distribution
    event_time_new[0] = rgeom1(m);
    
    // draw remaining loci from Markov process
    for (int l=1; l<L; ++l) {
      
      // copy values from previous locus
      event_time_new[l] = event_time_new[l-1];
      migrate_target_new[l] = migrate_target_new[l-1];
      
      // draw recombination events, leading to new draw from migration times
      double prob_recom = 1 - exp(-rho*(event_time_new[l] + 1)*delta[l]);
      recom_new[l] = rbernoulli1(prob_recom);
      if (recom_new[l]) {
        int split_time = sample2(0, event_time_new[l]);
        recom_time_new[l] = split_time;
        event_time_new[l] = split_time + rgeom1(m);
        migrate_target_new[l] = next_migrant++;
      }
    }
    
  }
  
  // return list
  return Rcpp::List::create(Rcpp::Named("recom_new") = recom_new,
                            Rcpp::Named("recom_time_new") = recom_time_new,
                            Rcpp::Named("event_time_new") = event_time_new,
                            Rcpp::Named("coalescence_new") = coalescence_new,
                            Rcpp::Named("coalesce_target_new") = coalesce_target_new,
                            Rcpp::Named("migrate_target_new") = migrate_target_new);
  */
  
  return Rcpp::List::create(Rcpp::Named("foo") = -9);
}
