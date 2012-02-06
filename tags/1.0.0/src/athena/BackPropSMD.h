// BackPropSMD.h
#ifndef __BACKPROPSMD_H__
#define __BACKPROPSMD_H__

#include "Terminals.h"
#include <Dataset.h>
#include "BackPropTree.h"
#include "BackProp.h"

///
/// Runs back propagation on the usual GENN network representation.
/// Makes a few assumptions on the how back propagation will work in 
/// this case.  
///
class BackPropSMD:public BackProp{

  public:
  
    BackPropSMD();
    
    /// Receives a postfix stack representation of the neural network
    int runBackProp(vector<TerminalSymbol*> & postfix_stack, Dataset * set);
    
    /// Sets whether to use stochastic
    void useBatch(bool tf){runBatch = tf;}

  private:
    void initialize();
    float runBatchEpoch(Dataset* set);
    float runStochasticEpoch(Dataset* set);
    double calculateMSE(Dataset* set);

    bool runBatch;
    int num_epochs, stochastic_size;

};

#endif
