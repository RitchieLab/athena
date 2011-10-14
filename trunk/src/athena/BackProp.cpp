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
