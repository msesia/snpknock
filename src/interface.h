#ifndef INTERFACE_H
#define INTERFACE_H

#include <RcppArmadillo.h>
#include <vector>
#include <iostream>
#include "dmc_old.h"
#include "dmc.h"
#include "hmm_old.h"
#include "hmm.h"
#include "haplotypes.h"
#include "genotypes.h"

typedef std::vector< double > vector;
typedef std::vector< std::vector<double> > vector2;
typedef std::vector< vector2 > vector3 ;
typedef std::vector< int > ivector;
typedef std::vector< std::vector<int> > ivector2;

vector  numToVec(const Rcpp::NumericVector &);
ivector numToIntVec(const Rcpp::IntegerVector & v);
vector2 numToVec2(const Rcpp::NumericVector &, int);
vector3 numToVec3(const Rcpp::NumericVector &, int, int);
ivector2 numToIntVec2(const Rcpp::IntegerVector &, int);

#endif
