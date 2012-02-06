//BackProp.h

#ifndef __BACKPROP_H__
#define __BACKPROP_H__

#include "Terminals.h"
#include "BackPropTree.h"
#include <Dataset.h>

class BackProp{

  public:
  
    BackProp();
  
    virtual ~BackProp(){}
  
    /// Receives a postfix stack representation of the neural network
    virtual int runBackProp(vector<TerminalSymbol*> & postfix_stack, Dataset * set)=0;
  
    /// Sets learning rate for the network
    void setLearningRate(double new_beta){learnRate = new_beta;}
  
    /// Returns learning rate
    double getLearningRate(){return learnRate;}
    
    /// Sets the number of iterations for running back propagation
    void setNumEpochs(int n_epochs){num_epochs = n_epochs;}
    
    /// Returns the number of iterations
    int getNumEpochs(){return num_epochs;}
    
    /// Sets threshold for stopping back propagation early
    void setBPThreshold(double new_thresh){threshold = new_thresh;}
    
    /// Returns the threshold
    double getBPThreshold(){return threshold;}
    
    /// Returns the list of weights in network
    virtual vector<float> getWeights(){return bptree.getWeights();}    

  protected:
    void initialize();
    
    virtual double calculateMSE(Dataset* set);
    
    BackPropTree bptree;
    
    double threshold, learnRate;
    int num_epochs, stochastic_size;

};

#endif
