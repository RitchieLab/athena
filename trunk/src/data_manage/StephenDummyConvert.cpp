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
#include "StephenDummyConvert.h"

#include <iostream>
using namespace std;

namespace data_manage
{

StephenDummyConvert::StephenDummyConvert()
{
}

StephenDummyConvert::~StephenDummyConvert()
{
}


///
/// Converts genotypes to ott dummy ones
///
void StephenDummyConvert::convert_genotypes(Dataholder* holder){

  Individual* ind;

  unsigned int curr_loc;
  unsigned int num_loci = holder->num_genos();

  for(unsigned int curr_ind = 0; curr_ind < holder->num_inds(); curr_ind++){

    ind = holder->get_ind(curr_ind);

    for(curr_loc=0; curr_loc < num_loci; curr_loc++){
      if(ind->get_genotype(curr_loc) < 3){
        ind->set_genotype(curr_loc, ind->get_genotype(curr_loc)-1);
      }

    }
  }

  // when it is ott dummy encoded there are 2 variables for each genotype
  // in this case we still only use one variable
  holder->ott_dummy_encoding(false);
}


}
