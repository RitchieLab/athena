// BackPropAnnie.h
#ifndef __BACKPROPANNIE_H__
#define __BACKPROPANNIE_H__

#include "Terminals.h"
#include <Dataset.h>
#include "BackPropTree.h"
#include "BackProp.h"
#include "BackPropAnnieTree.h"

///
/// Runs back propagation on the usual GENN network representation.
/// Makes a few assumptions on the how back propagation will work in 
/// this case.  
///

class BackPropAnnie:public BackProp{

  public:
  
    BackPropAnnie();
    
    /// Receives a postfix stack representation of the neural network
    int runBackProp(vector<TerminalSymbol*> & postfix_stack, Dataset * set);

    /// Returns the list of weights in network
    vector<float> getWeights(){return bpAnnieTree.getWeights();}
    
    /// Returns optimized score
    float getOptimizedScore(){return optimized_score;}

  private:
    
    void initialize();
    
    BackPropAnnieTree bpAnnieTree;

    double calculateMSE(Dataset* set);
    float optimized_score;
};

#endif
