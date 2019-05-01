/*
  This file is part of SNPknock.

    Copyright (C) 2017 Matteo Sesia

    SNPknock is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    SNPknock is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with SNPknock.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef O_KNOCKOFF_DMC_H
#define O_KNOCKOFF_DMC_H

#include <vector>
#include <random>
#include "utils.h"
#include <iostream>

typedef std::vector< std::vector<double> > matrix;

class oKnockoffDMC {
 public:
  oKnockoffDMC(const std::vector<double> & _initP, const std::vector< matrix > & _Q, int seed);
  ~oKnockoffDMC();
  std::vector<int> sample(const std::vector<int> & X);
  std::vector< std::vector<int> > sample(const std::vector<std::vector<int> > & X);
 private:
  std::vector<double> initP;
  std::vector< matrix > Q;
  int p, K;
  double tempQ;
  std::vector<double> Z, Z_old, W;
  std::vector<int> Xt;
  std::mt19937 gen;
  std::random_device rd;
  std::uniform_real_distribution<> dis; 
};

#endif
