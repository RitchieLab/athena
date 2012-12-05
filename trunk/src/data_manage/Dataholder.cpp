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
#include "Dataholder.h"
#include "Stringmanip.h"

namespace data_manage
{


Dataholder::Dataholder()
{
  any_missing = false;
  ott_encoded = false;
  max_locus = 0;
  split_num = -1;
}


///
/// Destructor frees all memory
///
Dataholder::~Dataholder()
{
  vector<Individual*>::iterator iter;
  for(iter = inds.begin(); iter != inds.end(); iter++)
    delete *iter;
}



///
/// Adds individual to set
/// @param ind Individual to add
///
void Dataholder::add_ind(Individual& ind){
  Individual* new_ind = new Individual(ind);
  inds.push_back(new_ind);
  inds_map[new_ind->get_id()] = new_ind;
}


///
/// Adds individual to set
/// @param ind Individual to add
///
void Dataholder::add_ind(Individual* ind){
  inds.push_back(ind);
  inds_map[ind->get_id()] = ind;
}

///
/// Adds default snp names to set
///
void Dataholder::add_default_snps(){
   unsigned int total = num_genos();
   for(unsigned int i=1; i<=total; i++){
     add_geno_name(Stringmanip::itos(i));
   }
}

///
/// Adds default covariate names to holder
///
void Dataholder::add_default_covars(){
    unsigned int total = inds[0]->num_covariates();
    for(unsigned int i=1; i<=total; i++){
      add_covar_name(Stringmanip::itos(i));  
    }
}


}
