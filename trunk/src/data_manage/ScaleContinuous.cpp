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
#include "ScaleContinuous.h"

#include <sstream>
#include <iostream>
using namespace std;

namespace data_manage
{

ScaleContinuous::ScaleContinuous()
{
}

ScaleContinuous::~ScaleContinuous()
{
}


///
/// For the continuous variable in question, the values are all
/// scaled by dividing by the largest value in the variable
/// @param holder Dataholder with all dat
/// @param var_index Variable to be scaled
///
void ScaleContinuous::adjust_contin(Dataholder* holder, unsigned int var_index){


  unsigned int curr_ind;
  Individual* ind;

  float max_value = holder->get_ind(0)->get_status();
  float covarmin = holder->get_ind(0)->get_status();
  float covaradjust = 0;
  
  for(curr_ind = 0; curr_ind < holder->num_inds(); curr_ind++){
    if(holder->get_ind(curr_ind)->get_covariate(var_index) == holder->get_missing_covalue()){
      continue;
    }
    if(holder->get_ind(curr_ind)->get_covariate(var_index) > max_value)
      max_value = holder->get_ind(curr_ind)->get_covariate(var_index);
    if(holder->get_ind(curr_ind)->get_covariate(var_index) < covarmin)
      covarmin = holder->get_ind(curr_ind)->get_covariate(var_index);
  }


  // when minimum is positive number use statmin as zero
  if(covarmin < 0){
    covaradjust = -covarmin;
    max_value = max_value + covaradjust;
  }

  // divide all values by largest and set in dataholder
  for(curr_ind=0; curr_ind < holder->num_inds(); curr_ind++){
    if(holder->get_ind(curr_ind)->get_covariate(var_index) == holder->get_missing_covalue())
      continue;
    ind = holder->get_ind(curr_ind);
    ind->set_covariate(var_index, (ind->get_covariate(var_index)+covaradjust)/max_value);
  }

}


///
/// Divides all status values by the maximum value so they
/// will be scaled from zero to one.
/// @param holder Dataholder with all data
///
void ScaleContinuous::adjust_status(Dataholder* holder){
  statmax = holder->get_ind(0)->get_status();
  statmin = holder->get_ind(0)->get_status();
  statadjust = 0;
  
  unsigned int curr_ind;
  Individual* ind;
  
  float status;
  
  for(curr_ind=0; curr_ind < holder->num_inds(); curr_ind++){
    status = holder->get_ind(curr_ind)->get_status();
    if(status > statmax)
      statmax = status;
    if(status < statmin)
      statmin = status;
  }
  
  // when minimum is positive number use statmin as zero
  if(statmin < 0){
    statadjust = -statmin;
    statmax = statmax + statadjust;
  }
  
  // divide all values by max and set status to that
  for(curr_ind=0; curr_ind < holder->num_inds(); curr_ind++){
    ind=holder->get_ind(curr_ind);
    ind->set_status((ind->get_status()+statadjust)/statmax);
  }  
}


///
/// Returns string that gives information on scaling performed
/// @return string
///
string ScaleContinuous::output_scale_info(){
  
  stringstream ss;
  
  ss << "ScaleMax=" << statmax << " StatusAdjust=" << statadjust << std::endl;
  return ss.str();
}

}
