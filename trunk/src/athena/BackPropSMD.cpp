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
// BackPropSMD.cpp

#include "BackPropSMD.h"
#include <RandomNoReplace.h>

///
/// Constructor
///
BackPropSMD::BackPropSMD(){
  initialize();
}

///
/// Initializer
///
void BackPropSMD::initialize(){
  BackProp::initialize();
  // size of a stochastic epoch
  stochastic_size = 10;
  // whether it should be run in batch or stochastic format
  runBatch = true;
}

///
/// Runs all individuals through as an epoch
/// @param set Dataset containing individuals
/// @return mean squared error for the set
/// 
float BackPropSMD::runBatchEpoch(Dataset* set){
  for(unsigned int indindex=0; indindex < set->num_inds(); indindex++){
    bptree.backPropagate((*set)[indindex]);
  }
 
  return calculateMSE(set);
}


///
/// Runs a sample of the individuals in the set
/// rather than the entire set.  If the learning rate is set
/// appropriately this can be as effective as running the batch and 
/// quicker.  This should be the case in cases where there
/// is duplicate data in the set.
/// @param 
///
float BackPropSMD::runStochasticEpoch(Dataset* set){
  
  vector<int> samples(stochastic_size, 0);
  
  RandomNoReplace::SampleWithoutReplacement(set->num_inds(), stochastic_size,
    samples);
  
  for(int i=0; i<stochastic_size; i++){
    bptree.backPropagate((*set)[samples[i]]);
  }
  
  return calculateMSE(set);
}

///
/// Receives a postfix stack representation of the neural network.
/// Converts to a tree representation.  Runs back propagation and changes the
/// weights in the stack to match the result of back propagation.
/// @param postfix_stack vector of TerminalSymbols representing neural network
/// @param set Dataset used for running backpropagation
/// @returns number of epochs
///
int BackPropSMD::runBackProp(vector<TerminalSymbol*> & postfix_stack, Dataset* set){

  // create the back propagation tree that will be used to adjust the weights 
  // in the original neural network
  bptree.createTree(postfix_stack);
  
  float mse, lastmse=0.0;
  int numEpochs = 0;
  
  // will run a number of epochs (the epochs will either be everyone in the set)
  // or a subset of the total
  if(runBatch){
    for(int epoch=0; epoch < num_epochs; epoch++){
      numEpochs++;
      mse = runBatchEpoch(set);
      if(fabs(mse-lastmse) < 0.00001)
	break;
      else
        lastmse = mse;
    }
  }
  else{
    
    for(int epoch=0; epoch < num_epochs; epoch++){
      mse = runStochasticEpoch(set);
      if(fabs(mse-lastmse) < 0.0001)
	break;
      else
        lastmse = mse;
    }
  }

  return numEpochs;

}

///
/// Calculates mean-squared error on current backprop tree using 
/// current weights.
/// @param set Dataset used for running back prop
/// @return mse
///
double BackPropSMD::calculateMSE(Dataset* set){
  
  double total=0.0, output;
  for(unsigned int indindex=0; indindex < set->num_inds(); indindex++){
    output =  bptree.feedForward((*set)[indindex]);
    total += (output - (*set)[indindex]->status()) * (output - (*set)[indindex]->status());
  }
  return  total/set->num_inds();
  
}


