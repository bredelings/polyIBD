
#pragma once

#include <Rcpp.h>

//------------------------------------------------
// define emmission probability lookup table
std::vector< std::vector< std::vector< std::vector< std::vector<double> > > > > define_emmission_lookup(int m_max, int L, std::vector<double> &p_vec);

