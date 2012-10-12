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
//BackProp.cpp
#include "BackProp.h"


///
/// Constructor
///
BackProp::BackProp(){
  initialize();
}


///
/// Initialize parameters
///
void BackProp::initialize(){
    // learning rate
  learnRate = 0.3;
  // value of target mean-square error (will stop when that value is reached)
  threshold = 0.000001;
  // base number of epochs 
  num_epochs = 500;
}


///
/// Calculates mean-squared error on current backprop tree using 
/// current weights.
/// @param set Dataset used for running back prop
/// @return mse
///
double BackProp::calculateMSE(Dataset* set){
  
  double total=0.0, output;
  for(unsigned int indindex=0; indindex < set->num_inds(); indindex++){
    output =  bptree.feedForward((*set)[indindex]);
    total += (output - (*set)[indindex]->status()) * (output - (*set)[indindex]->status());
  }
  return  total/set->num_inds();
  
}
