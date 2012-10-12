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

#ifndef INDIVIDUAL_H_
#define INDIVIDUAL_H_

#include <vector>
#include <string>
#include <map>

namespace data_manage
{

/// Contains data and status for individuals in dataset
class Individual
{
public:
  Individual();

  ~Individual();

  /// status of individual
  inline float status(){return ind_status;}

  /// returns genotype or environmental value
  inline float value(unsigned int index){
    if(index<num_covars) return covariates[index];
    else return float(genotype[index+num_covars]);
  }

  inline int get_genotype(unsigned int index){return genotype[index];}
  
  inline void set_genotype(unsigned int index, char value){genotype[index]=value;}

  inline void set_status(float stat){ind_status = stat;}

  inline float get_status(){return ind_status;}

  /// appends genotype to genotype lsit
  inline void add_genotype(int geno){genotype.push_back(geno);num_loci++;}

  /// appends environmental variable to list
  void add_covariate(float val){covariates.push_back(val);num_covars++;}

  /// returns covariate value
  float get_covariate(unsigned int index){return covariates[index];}

  /// sets covariate value
  void set_covariate(unsigned int index, float val){covariates[index] = val;}

  /// returns number of genotypes
  inline unsigned int num_genotypes(){return num_loci;}

  /// returns number of environmental variables
  inline unsigned int num_covariates(){return num_covars;}

  /// sets genotype vector to be equal to vector passed
  inline void set_all_genotypes(std::vector<char>& genos){genotype = genos; num_loci = genos.size();}

  /// returns ID value as string
  std::string get_id(){return id;}
  
  /// sets ID for individual 
  void set_id(std::string identification){id=identification;}
  
private:
  float ind_status;
  unsigned int num_loci, num_covars;
  std::string id;
  std::vector<float> covariates;
  std::vector<char> genotype;

};

}

#endif /*INDIVIDUAL_H_*/

