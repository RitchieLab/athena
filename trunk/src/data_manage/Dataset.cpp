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

#include "Dataset.h"

using namespace std;

namespace data_manage
{

Dataset::Dataset()
{
  sstotal = 0;
}

Dataset::~Dataset()
{
}


///
/// Appends pointer list to end of current list
/// @param new_inds vector of Individual pointers
///
void Dataset::add_inds(vector<Individual* >& new_inds){
  inds.insert(inds.end(), new_inds.begin(), new_inds.end());
}


///
/// Outputs set displaying the status and genotypes for use in verifying
/// the sets
/// @param os ostream to write to
/// @param d Dataset containes current individuals
/// @return ostream
///
ostream& operator<<(ostream& os, Dataset& d){
  Individual * ind;
  for(unsigned int i=0; i<d.num_inds(); i++){
    ind = d[i];
    os << ind->get_status();
    for(unsigned int loc=0; loc < d.num_genos(); loc++){
      os << " " << ind->get_genotype(loc);
    }
    // output continuous after genotypes
    for(unsigned int cov=0; cov < d.num_covariates(); cov++){
      os << " " << ind->get_covariate(cov);
    }
    
    os << endl;
  }
  
  return os;
}


///
/// Calculates SSTotal on this set
///
void Dataset::calc_sstotal(){
  
  float diff=0.0;
  float meanval;
  float stat_total = 0.0;
  
  for(unsigned int i=0; i < inds.size(); i++){
    stat_total += inds[i]->get_status();
  }
  
  meanval = stat_total/inds.size();
  
  for(unsigned int i=0; i<inds.size(); i++){
    diff = diff + (inds[i]->get_status() - meanval) * (inds[i]->get_status() - meanval);
  }
  
  sstotal = diff;
}


}
