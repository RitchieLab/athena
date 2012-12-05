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
//BackPropAnnieTree.h

#ifndef __BACKPROPANNIETREE_H__
#define __BACKPROPANNIETREE_H__


#include "Terminals.h"
#include "ExpressionTree.h"
#include "Tree.hh"
#include <Dataset.h>
#include <annie.h>

///
/// Tree based representation that works well with back propagation.
///

class BackPropAnnieTree{

  public:
    BackPropAnnieTree();
    BackPropAnnieTree(vector<TerminalSymbol*>& postfix_stack, data_manage::Dataset* set);
    ~BackPropAnnieTree();
  
    void createTree(vector<TerminalSymbol*>& postfix_stack, data_manage::Dataset* set);
    float feedForward(Individual* ind);
    
    /// Returns weights from current network
    vector<float> getWeights();
  
    /// Train network using annie library
    void trainNetwork(Dataset* set, int epochs);
  
  private:
    
    void initialize();
    
    struct Neuron{
      Neuron(){output=0.0;delta=0.0;prevWeight=0.0;currWeight=0.0;
        constantWeightValue=NULL;adjusted_depth=0;}

      inline void setAnnieLink(int layer, int neuron, int link){
        annie_layer = layer;
        annie_neuron = neuron;
        annie_link = link;
      }
      
      TerminalSymbol* neuronTerminal;
      // original constant terminal value for this neuron
      // (going up to next one) -- will be blank on the output neuron
      TerminalSymbol* constantWeightValue;
      
      // annie network stores connections with network above in the network
      // unlike this implementation
      int annie_layer; // layer for connection to in annie network
      int annie_neuron; // index of the neuron to connect to 
      int annie_link; // index of neuron connecting to it
      
      double prevWeight;
      double currWeight;
      double output;
      double delta;
      int adjusted_depth;
    };
    
    struct inputValue{
      int index;
      bool isGenotype;
    };
    
    tree<Neuron> bptree;
    float beta;
    int maxDepth;
    
    void addNeuron(ExpressTreeIter& exIter, tree<Neuron>::iterator& bpIter,
      ExpressionTree& extree, tree<Neuron>& bptree);
    
    void setMultiNetworkConnections(tree<Neuron>::iterator& bpIter,
      tree<Neuron>& bptree, vector<int>& neuronIndexes, vector<int>& destIndexes,
      annie::MultiLayerNetwork* net);
    
    void constructTrainingSet(data_manage::Dataset* set);
    
    annie::MultiLayerNetwork* net;
    annie::TrainingSet* trainSet;
    vector<inputValue> inputIndexes;    
    

};


#endif
