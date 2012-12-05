/*
Copyright Marylyn Ritchie 2011

This file is part of ATHENA.

ATHENA is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

ATHENA is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with ATHENA.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef DATASET_H_
#define DATASET_H_

#include "Individual.h"
#include <iostream>

namespace data_manage
{

///
///	Dataset contains pointers to individual objects.  It can represent the entire
/// data set or can be used for a permuted set or a training, testing or crossvalidation
/// interval.
///
class Dataset
{

  friend std::ostream& operator<<(std::ostream& os, Dataset& d);

public:
  Dataset();
  Dataset(float comissing, unsigned int genomissing){
      missing_covalue = comissing;
      missing_genotype = genomissing;
      sstotal = 0;
  }
  ~Dataset();

  inline unsigned int num_inds(){return inds.size();}

  inline void add_ind(Individual* ind){inds.push_back(ind);}

  inline Individual* operator[](unsigned int index){return inds[index];}

  inline Individual* get_ind(unsigned int index){return inds[index];}
  
  /// Returns number of genotypes in set
  inline unsigned int num_genos(){return inds[0]->num_genotypes();}
  
  /// Returns number of covariates in set
  inline unsigned int num_covariates(){return inds[0]->num_covariates();}

  void add_inds(std::vector<Individual* >& new_inds);

  inline float get_missing_covalue(){return missing_covalue;}
  
  inline void set_missing_covalue(float miss){missing_covalue = miss;}
  
  inline int get_missing_genotype(){return missing_genotype;}
  
  inline void set_missing_genotype(int miss){missing_genotype = miss;}
  
  /// returns SStotal for this set (used in r-squared calculations)
  float get_sstotal(){return sstotal;}
  
  /// Calculates SStotal
  void calc_sstotal();
  
private:
  std::vector<Individual*> inds;
  float missing_covalue, sstotal;
  int missing_genotype;

};

}

#endif /*DATASET_H_*/
