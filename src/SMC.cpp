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
  // population that the genetic element in sample i at locus j occupies
  vector<vector<int>> x(n);
  for (int i=0; i<n; ++i) {
    x[i] = vector<int>(L, i);
  }
  
  // store whether lineage has coalesced, which lineage it coalesced with, and
  // when coalescence occured
  vector<vector<bool>> coalesced(n, vector<bool>(L, false));
  vector<vector<int>> coalesced_with(n, vector<int>(L));
  vector<vector<int>> coalesced_time(n, vector<int>(L));
  
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
        
        // skip if already coalesced
        if (coalesced[i][l]) {
          continue;
        }
        
        // if first time seeing the value x[i][l] as we scan across loci then we
        // are essentially jumping to a new host, therefore draw new parent. If
        // we have seen x[i][l] before then we have already drawn a new parent
        // for the genetic material in this host. which will be stored in the
        // recom_map. Therefore, draw new or existing parent based on
        // recombination probability and distance since previous value
        int new_value = 0;
        if (recom_map.count(x[i][l]) != 1) {  // first time seeing this value
          
          // draw new parent
          new_value = sample2(0,N-1);
          
        } else {  // seen before
          
          // get probability of recombination
          int delta = loci[l] - recom_map[x[i][l]].first;  // distance since previous value
          double prob_break = 1 - exp(-rho*delta);
          
          // if recombination then draw new parent, otherwise same parent as
          // previous locus with this value
          if (rbernoulli1(prob_break)) {
            new_value = sample2(0,N-1);
          } else {
            new_value = recom_map[x[i][l]].second;  // the stored parent
          }
          
        }
        
        // update recom_map and x
        recom_map[x[i][l]] = {loci[l], new_value};
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
            coalesced_with[i][l] = j;
            coalesced_time[i][l] = t;
            if (j == 0) {
              coalesced[j][l] = true;
              coalesced_with[j][l] = -1;
              coalesced_time[j][l] = t;
            }
            
          }
        }
      }
    }
    
  } // end loop through time
  
  // return list
  return Rcpp::List::create(Rcpp::Named("coalesced") = coalesced,
                            Rcpp::Named("coalesced_with") = coalesced_with,
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
      
      // draw timing of recombination event uniformly over branch lengths
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

//------------------------------------------------
// draw from SMC under wright-fisher model by looping through discrete
// generations
// [[Rcpp::export]]
Rcpp::List wrightfisher_SMC_naive_cpp(Rcpp::List args) {
  
  // get inputs
  int n = Rcpp_to_int(args["n"]);
  vector<int> loci = Rcpp_to_vector_int(args["loci"]);
  int L = int(loci.size());
  int N = Rcpp_to_int(args["N"]);
  double rho = Rcpp_to_double(args["rho"]);
  int generations = Rcpp_to_int(args["generations"]);
  
  // get distance between loci
  vector<int> delta(L);
  for (int l=1; l<L; ++l) {
    delta[l] = loci[l] - loci[l-1];
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
  vector<vector<int>> coalesced_time(n, vector<int>(L));
  
  // loop backwards through time
  for (int t=0; t<generations; ++t) {
    
    // implement recombination/independent sampling
    for (int i=0; i<n; ++i) {
      int prev_locus = x[i][0];
      int prev_locus_new_value = 0;
      for (int l=0; l<L; ++l) {
        
        // skip if already coalesced
        if (coalesced[i][l]) {
          continue;
        }
        
        // draw recombination
        double prob_recom = 1;
        if (l > 0 && x[i][l] == prev_locus) {
          prob_recom = 1 - exp(-rho*delta[l]);
        }
        if (rbernoulli1(prob_recom)) {
          prev_locus_new_value = sample2(0,N-1);
          prev_locus = x[i][l];
        }
        x[i][l] = prev_locus_new_value;
      }
    }
    
    // check for coalescence
    for (int l=0; l<L; ++l) {
      for (int i=1; i<n; ++i) {
        
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
            
          }
        }
      }
    }
    
  }  // end loop over time
  
  // return list
  return Rcpp::List::create(Rcpp::Named("coalesced") = coalesced,
                            Rcpp::Named("coalesced_time") = coalesced_time);
}

//------------------------------------------------
// draw from conditional SMC under wright-fisher model
// [[Rcpp::export]]
Rcpp::List wrightfisher_SMC_conditional_cpp(Rcpp::List args) {
  
  // get inputs
  int n = Rcpp_to_int(args["n"]);
  vector<int> loci = Rcpp_to_vector_int(args["loci"]);
  int L = int(loci.size());
  int N = Rcpp_to_int(args["N"]);
  double rho = Rcpp_to_double(args["rho"]);
  
  // get distance between loci
  vector<int> delta(L);
  for (int l=1; l<L; ++l) {
    delta[l] = loci[l] - loci[l-1];
  }
  
  // cw = coalesce with. Integer value giving which other lineage a given
  // lineage coalesces to
  // ct = coalesce time. Time at which coalescence occurs
  vector<vector<int>> cw(L, vector<int>(n-1,-1));
  vector<vector<int>> ct(L, vector<int>(n-1,-1));
  
  // loop through samples
  vector<int> ct_sort;
  for (int i=0; i<(n-1); ++i) {
    
    // draw first locus from standard coalescent process
    ct_sort = ct[0];
    sort(ct_sort.begin(), ct_sort.begin()+i);
    int j = 0;
    int k = i+1;
    int coalesce_time = rgeom1(k/double(N));
    while (coalesce_time >= ct_sort[j] && k > 1) {
      j++;
      k--;
      while (ct_sort[j] == ct_sort[j-1] && k > 1) {
        j++;
        k--;
      }
      coalesce_time = ct_sort[j-1] + rgeom1(k/double(N));
    }
    ct[0][i] = coalesce_time;
    
    // choose random lineage to coalesce with
    int coalesce_with = 1;
    int tmp1 = sample2(0,k-1);
    for (int j=0; j<i; ++j) {
      if (ct[0][j] > coalesce_time) {
        if (tmp1 == 0) {
          coalesce_with = j+2;
          break;
        }
        tmp1--;
      }
    }
    cw[0][i] = coalesce_with;
    
    // loop through remaining loci
    for (int l=1; l<L; ++l) {
      
      // get maximum coalescent time ignoring lineage i
      int t_max = 0;
      for (int j=0; j<i; ++j) {
        if (ct[l][j] > t_max) {
          t_max = ct[l][j];
        }
      }
      
      // use t_max to calculate overhang
      int overhang = (ct[l-1][i] > t_max) ? ct[l-1][i] - t_max : 0;
      
      // draw whether there is recombination event since previous locus
      double prob_recom = 1 - exp(-rho*(ct[l-1][i] + overhang + 1)*delta[l]);
      if (rbernoulli1(prob_recom)) {
        
        // draw recombination split time uniformly over overhang or this branch
        double prob_overhang = 0;
        if (overhang > 0) {
          prob_overhang = overhang/double(overhang + ct[l-1][i]);
        }
        int split_time = 0;
        if (rbernoulli1(prob_overhang)) {
          split_time = sample2(t_max, ct[l-1][i]);
        } else {
          split_time = sample2(0, ct[l-1][i]);
        }
        
        // get number of extant lineages k at split_time
        ct_sort = ct[l-1];
        sort(ct_sort.begin(), ct_sort.begin()+i);
        int k = i+1;
        for (int j=0; j<i; ++j) {
          if (ct_sort[j] < split_time) {
            k--;
          }
        }
        
        // draw time at which lineage rejoins tree
        coalesce_time = split_time + rgeom1(k/double(N));
        int j = i+1-k;
        while (coalesce_time >= ct_sort[j] && k > 1) {
          j++;
          k--;
          while (ct_sort[j] == ct_sort[j-1] && k > 1) {
            j++;
            k--;
          }
          coalesce_time = ct_sort[j-1] + rgeom1(k/double(N));
        }
        
        // choose random lineage to coalesce with
        int tmp1 = sample2(0,k-1);
        coalesce_with = 1;
        for (int j=0; j<i; ++j) {
          if (ct[l][j] > coalesce_time) {
            if (tmp1 == 0) {
              coalesce_with = j+2;
              break;
            }
            tmp1--;
          }
        }
        
      }  // end bernoulli prob_recom
      
      // update coalescence
      cw[l][i] = coalesce_with;
      ct[l][i] = coalesce_time;
      
    }  // end loop over loci
  } // end loop over samples
  
  // return list
  return Rcpp::List::create(Rcpp::Named("coalesced_with") = cw,
                            Rcpp::Named("coalesced_time") = ct);
}

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
  int n = Rcpp_to_int(args["n"]);
  vector<int> loci = Rcpp_to_vector_int(args["loci"]);
  int L = int(loci.size());
  int N = Rcpp_to_int(args["N"]);
  double rho = Rcpp_to_double(args["rho"]);
  double mu = Rcpp_to_double(args["mu"]);
  
  // get distance between loci
  vector<int> delta(L);
  for (int l=1; l<L; ++l) {
    delta[l] = loci[l] - loci[l-1];
  }
  
  int next_migrant = 0;
  
  vector<vector<int>> event_time(L, vector<int>(n,-1));
  vector<vector<bool>> coalesced(L, vector<bool>(n,false));
  vector<vector<int>> coalesced_with(L, vector<int>(n,-1));
  vector<vector<int>> migrated_with(L, vector<int>(n,-1));
  
  // loop through samples
  vector<int> t_sort;
  for (int i=0; i<n; ++i) {
    
    // --- draw first locus from standard coalescent process ---
    
    // draw time at which lineage rejoins tree
    t_sort = event_time[0];
    sort(t_sort.begin(), t_sort.begin()+i);
    int new_time = 0;
    int k = i;
    for (int j=0; j<(i+1); ++j) {
      k = i - j;
      new_time += rgeom1(k/double(N) + mu);
      if (new_time < t_sort[j] || k == 0) {
        event_time[0][i] = new_time;
        double prob_coalesce = k/(k + N*mu);
        coalesced[0][i] = rbernoulli1(prob_coalesce);
        break;
      }
    }
    
    // choose random lineage to coalesce with
    if (coalesced[0][i]) {
      int tmp1 = sample2(0,k-1);
      for (int j=0; j<i; ++j) {
        if (event_time[0][j] > new_time) {
          if (tmp1 == 0) {
            coalesced_with[0][i] = j;
            break;
          }
          tmp1--;
        }
      }
    } else {
      migrated_with[0][i] = next_migrant++;
    }
    
    // --- draw remaining loci from Markov process ---
    
    // loop through remaining loci
    for (int l=1; l<L; ++l) {
      
      // copy over previous locus values
      event_time[l][i] = event_time[l-1][i];
      coalesced[l][i] = coalesced[l-1][i];
      coalesced_with[l][i] = coalesced_with[l-1][i];
      migrated_with[l][i] = migrated_with[l-1][i];
      
      // draw whether there is recombination event since previous locus
      double prob_recom = 1 - exp(-rho*(event_time[l-1][i] + 1)*delta[l]);
      if (rbernoulli1(prob_recom)) {
        
        // draw recombination split time uniformly over this branch
        int split_time = sample2(0, event_time[l-1][i]);
        
        // draw time at which lineage rejoins tree
        t_sort = event_time[l];
        sort(t_sort.begin(), t_sort.begin()+i);
        int k = i;
        int new_time = split_time;
        for (int j=0; j<(i+1); ++j) {
          if (t_sort[j] < split_time) {
            continue;
          }
          k = i - j;
          new_time += rgeom1(k/double(N) + mu);
          if (new_time < t_sort[j] || k == 0) {
            event_time[l][i] = new_time;
            double prob_coalesce = k/(k + N*mu);
            coalesced[l][i] = rbernoulli1(prob_coalesce);
            break;
          }
        }
        
        // choose random lineage to coalesce with
        if (coalesced[l][i]) {
          int tmp1 = sample2(0,k-1);
          for (int j=0; j<i; ++j) {
            if (event_time[l][j] > new_time) {
              if (tmp1 == 0) {
                coalesced_with[l][i] = j;
                break;
              }
              tmp1--;
            }
          }
        } else {
          migrated_with[l][i] = next_migrant++;
        }
        
      }  // end recom prob
      
    }  // end loop over loci
    
  } // end loop over samples
  
  // return list
  return Rcpp::List::create(Rcpp::Named("event_time") = event_time,
                            Rcpp::Named("coalesced") = coalesced,
                            Rcpp::Named("coalesced_with") = coalesced_with,
                            Rcpp::Named("migrated_with") = migrated_with);
  
}

//------------------------------------------------
// draw from conditional SMC under continent-island model
// [[Rcpp::export]]
Rcpp::List continentisland_SMC_conditional2_cpp(Rcpp::List args) {
  
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
