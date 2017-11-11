
#pragma once

#include <Rcpp.h>

//------------------------------------------------
// define emmission probability lookup table
std::vector< std::vector< std::vector< std::vector< std::vector<double> > > > > define_emmission_lookup(int m_max, int L, std::vector<double> &p_vec, const double e1, const double e2);

//------------------------------------------------
// update transition probability lookup table
void update_transition_lookup(std::vector< std::vector< std::vector<double> > > &transition_lookup, double f, double rho, int z_max);

//------------------------------------------------
// text
double forward_alg(std::vector< std::vector<double> > &frwrd_mat, const std::vector<int> &x, const std::vector< std::vector< std::vector<double> > > &emmission_lookup, const std::vector< std::vector< std::vector<double> > > &transition_lookup, double f);

//------------------------------------------------
// text
void backward_alg(std::vector< std::vector<double> > &bkwrd_mat, const std::vector<int> &x, const std::vector< std::vector< std::vector<double> > > &emmission_lookup, const std::vector< std::vector< std::vector<double> > > &transition_lookup);

//------------------------------------------------
// text
void get_IBD(std::vector< std::vector<double> > &IBD_mat, const std::vector< std::vector<double> > &frwrd_mat, const std::vector< std::vector<double> > &bkwrd_mat);
