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

#include "NormalizeContinuous.h"
#include <Descriptive.h>

#include <iostream>
using namespace std;

namespace data_manage
{

NormalizeContinuous::NormalizeContinuous()
{
}

NormalizeContinuous::~NormalizeContinuous()
{
}


///
/// Rescales the continuous variables by taking absolute difference
/// from the mean and dividing by the standard deviation
/// @param holder Dataholder with all dat
/// @param var_index Variable to be scaled
///
void NormalizeContinuous::adjust_contin(Dataholder* holder, unsigned int var_index){
  stat::Descriptive calculator;

  vector<float> values(holder->num_inds(), 0);

  unsigned int curr_ind;
  for(curr_ind = 0; curr_ind < holder->num_inds(); curr_ind++){
    values[curr_ind] = holder->get_ind(curr_ind)->get_covariate(var_index);
  }

  float mean_val = calculator.mean(values);
  float std_dev = calculator.standard_dev(values);

  Individual* ind;
  // subtract each value from the mean and divide by the standard deviation
  for(curr_ind=0; curr_ind < holder->num_inds(); curr_ind++){
    ind = holder->get_ind(curr_ind);
    ind->set_covariate(var_index, ((ind->get_covariate(var_index) - mean_val) / std_dev));
  }

}



///
/// Rescales the continuous variables by taking absolute difference
/// from the mean and dividing by the standard deviation
/// @param holder Dataholder with all dat
///
void NormalizeContinuous::adjust_status(Dataholder* holder){
  stat::Descriptive calculator;
  
  vector<float> values(holder->num_inds(), 0);
  
  unsigned int curr_ind;
  for(curr_ind=0; curr_ind < holder->num_inds(); curr_ind++){
    values[curr_ind] = holder->get_ind(curr_ind)->get_status();
  }
  meanval = calculator.mean(values);
  stddev = calculator.standard_dev(values);

  Individual* ind;
  // subtract each value from the mean and divide by the standard deviation
  for(curr_ind=0; curr_ind < holder->num_inds(); curr_ind++){
    ind=holder->get_ind(curr_ind);
    ind->set_status((ind->get_status()-meanval)/stddev);
  }
  
}

}
