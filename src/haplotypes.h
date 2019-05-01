#ifndef HAPLOTYPES_H
#define HAPLOTYPES_H

/*
Knockoffs for phased haplotypes
*/

#include <vector>
#include <random>
#include "utils.h"
#include <iostream>

typedef std::vector< std::vector<double> > matrix;
typedef std::vector< std::vector<int> > imatrix;

class GroupHaplotypes {
 public:
  GroupHaplotypes(const std::vector<double> & r, const matrix & alpha, const matrix & _theta, 
                  const std::vector<int> _groups, int seed);
  imatrix sample(const imatrix & X);
  std::vector<int> sample(const std::vector<int> & X);
 private:
  void sampleViterbi(const std::vector<int> & X);
  void knockoffMC(const std::vector<int> & H);
  void emission(const std::vector<int> & Hk);
  matrix theta, a;
  std::vector<double> b;
  matrix beta;
  double beta_const;
  std::vector<double> weightsEmit, weights;
  int nStates, p, nGroups;
  std::vector<int> H, Hk, Xk, groups;
  imatrix elements;
  // Partition function for Markov chain knockoffs
  std::vector<double> Z, Z_old;
  // Random number generation
  std::random_device rd;
  std::mt19937 gen;
  std::uniform_real_distribution<> dis;
  std::mt19937 gen2;
};

#endif
