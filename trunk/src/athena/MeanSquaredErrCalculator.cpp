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
#include "MeanSquaredErrCalculator.h"

using namespace std;

MeanSquaredErrCalculator::MeanSquaredErrCalculator(){
    reset();
    name = "Mean Squared Error";
}

void MeanSquaredErrCalculator::reset(){
   total_inds_tested = 0;
   squared_error_total = 0.0;
   sstotal = 0.0;
   stat_total.clear();
}

///
/// Adds score to running total within object
/// @param score
///
void MeanSquaredErrCalculator::add_ind_score(float score, float stat){

    float difference = score - stat;
    squared_error_total += difference * difference;
    total_inds_tested++;
    sstotal += stat;
    stat_total.push_back(stat);
}

///
/// Calculates and returns sstotal for use in calculating
/// R-squared value
/// @return sstotal
///
float MeanSquaredErrCalculator::get_constant(){
  float mean = sstotal / total_inds_tested;
  float diff=0.0;
  for(vector<float>::iterator iter=stat_total.begin(); iter != stat_total.end();
    ++iter){
    diff = diff + (*iter-mean) * (*iter-mean);
  }  
  sstotal = diff;
  return sstotal;
}
