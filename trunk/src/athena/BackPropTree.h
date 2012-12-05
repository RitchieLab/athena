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
//BackPropTree.h

#ifndef __BACKPROPTREE_H__
#define __BACKPROPTREE_H__


#include "Terminals.h"
#include "ExpressionTree.h"
#include "Tree.hh"

///
/// Tree based representation that works well with back propagation.
///

class BackPropTree{

  public:
    BackPropTree();
    BackPropTree(vector<TerminalSymbol*>& postfix_stack);
  
    void createTree(vector<TerminalSymbol*>& postfix_stack);
    float feedForward(Individual* ind);
    
    /// Run backpropagation on this single individual
    float backPropagate(Individual* ind);
    
    /// Sets beta and runs backprop on this individual
    inline float backPropagate(Individual* ind, float b){
      setBeta(b);
      return backPropagate(ind);
    }
    
    /// Sets beta (learning rate) for the tree
    inline void setBeta(float b){beta = b;}
    
    /// Returns weights from current network
    vector<float> getWeights();
  
  private:
    
    void initialize();
    
    struct Neuron{
      Neuron(){output=0.0;delta=0.0;prevWeight=0.0;currWeight=0.0;
        constantWeightValue=NULL;}

      TerminalSymbol* neuronTerminal;
      // original constant terminal value for this neuron
      // (going up to next one) -- will be blank on the output neuron
      TerminalSymbol* constantWeightValue; 
      double prevWeight;
      double currWeight;
      double output;
      double delta;
    };
    
    tree<Neuron> bptree;
    float beta;
    
    void addNeuron(ExpressTreeIter& exIter, tree<Neuron>::iterator& bpIter,
      ExpressionTree& extree, tree<Neuron>& bptree);
};


#endif
