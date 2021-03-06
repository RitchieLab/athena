// BackPropAnnieAnnie.cpp

#include "BackPropAnnie.h"
#include "BackPropAnnieTree.h"
#include <RandomNoReplace.h>

///
/// Constructor
///
BackPropAnnie::BackPropAnnie(){
}

///
/// Receives a postfix stack representation of the neural network.
/// Converts to a tree representation.  Runs back propagation and changes the
/// weights in the stack to match the result of back propagation.
/// @param postfix_stack vector of TerminalSymbols representing neural network
/// @param set Dataset used for running backpropagation
/// @returns number of epochs
///
int BackPropAnnie::runBackProp(vector<TerminalSymbol*> & postfix_stack, Dataset* set){

  // create the back propagation tree that will be used to adjust the weights 
  // in the original neural network
  bpAnnieTree.createTree(postfix_stack, set);
  
  float mse, lastmse=1000.0;

  int numepochs=0;
  
  for(int i=0; i<100; i++){
    bpAnnieTree.trainNetwork(set, 1);
    numepochs += 1;
    bpAnnieTree.getWeights();
    mse = calculateMSE(set);
    // original difference was 0.0001
    // stop when Mean-squared error is not changing or the mean-squared error is getting bigger
    if(fabs(mse-lastmse)<threshold || mse > lastmse){
      break;
    }
    else{
      lastmse = mse;
    }
  }  
  
  optimized_score = mse;
  return numepochs;
}

///
/// Calculates mean-squared error on current backprop tree using 
/// current weights.
/// @param set Dataset used for running back prop
/// @return mse
///
double BackPropAnnie::calculateMSE(Dataset* set){
  double total=0.0, output;
  for(unsigned int indindex=0; indindex < set->num_inds(); indindex++){
    output =  bpAnnieTree.feedForward((*set)[indindex]);
    total += (output - (*set)[indindex]->status()) * (output - (*set)[indindex]->status());
  }
  return  total/set->num_inds();
  
}


