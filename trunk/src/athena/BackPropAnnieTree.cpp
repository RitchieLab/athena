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
#include "BackPropAnnieTree.h"
#include <deque>
#include <sstream>

using namespace annie;

BackPropAnnieTree::BackPropAnnieTree(){
  initialize();
}


BackPropAnnieTree::~BackPropAnnieTree(){
// cout << "trainSet in destructor=" << trainSet << endl;
  if(trainSet != NULL){
    delete trainSet;
  }
  if(net != NULL){
    delete net;
  }
}


///
/// Alternative constructor.  Creates a back prop tree from the 
/// postfix stack.
/// @param postfix_stack Vector containing terminals which will be converted
/// into a more back propagation friendly format.
/// @param set Dataset
///
BackPropAnnieTree::BackPropAnnieTree(vector<TerminalSymbol*>& postfix_stack, Dataset* set){
  initialize();
  createTree(postfix_stack, set);
}


///
/// Initializes object
///
void BackPropAnnieTree::initialize(){
  beta = 0.3;
  maxDepth = 0;
  net = NULL;
  trainSet = NULL;
}


///
/// Recursively adds nodes to the backprop tree from the expression tree.
/// @param exIter iterator pointing to current neuron in the expression tree
/// @param bpIter iterator pointing to current neuron in backprop tree
/// @param extree ExpressionTree 
/// @param bptree Backpropagation tree
///
void BackPropAnnieTree::addNeuron(ExpressTreeIter& exIter, tree<Neuron>::iterator& bpIter,
  ExpressionTree& extree, tree<Neuron>& bptree){

  deque<float> args;

  for(int child=0; child < extree.number_of_children(exIter); child++){
    Neuron newNode;
    
    ExpressTreeIter weightIter, childIter, exNodeIter;
    // each child will be a weight terminal symbol
    weightIter = extree.child(exIter, child);
    bool any_children = false;
   
    // now get the actual weight value and variable/neuron
    for(int wchild=0; wchild < extree.number_of_children(weightIter); wchild++){
      childIter=extree.child(weightIter, wchild);
      if(childIter->el->get_term_type() == TerminalSymbol::Constant){
        newNode.constantWeightValue = childIter->el;
        newNode.currWeight = childIter->el->evaluate(args);
      }
      // in this case it is either a variable or another neuron
      else{
        newNode.neuronTerminal = childIter->el;
        exNodeIter = childIter;
        if(extree.number_of_children(childIter) > 0)
          any_children=true;
      }
    }
    tree<Neuron>::iterator newBPnode = bptree.append_child(bpIter, newNode);
    
    // check max depth here 
    if(bptree.depth(newBPnode) > maxDepth){
      maxDepth = bptree.depth(newBPnode);
    }
    
    // if any children of this neuron (it is a PADD) then recursively
    // call this function to continue building backprop tree
    if(any_children){
      addNeuron(exNodeIter, newBPnode, extree, bptree);
    }
  }
  
}


///
/// Calculates feed-forward result on individual specified
/// @param ind Individual whose values will be 
/// @return The output of network
///
float BackPropAnnieTree::feedForward(Individual* ind){
  
  ContinVariable::set_ind(ind);
  GenotypeTerm::set_ind(ind);
  
  tree<Neuron>::iterator childIter;
  deque<float> args;
  
  // need to begin at bottom of tree and then progress up
  // record output at each neuron
  for(tree<Neuron>::post_order_iterator ffIter=bptree.begin_post();
     ffIter != bptree.end_post(); ffIter++){

    args.clear();
    
    // need to set output -- each neurons output is set by 
    // taking the outputs of neurons below -- adding together
    for(unsigned int child=0; child < bptree.number_of_children(ffIter); child++){
      childIter = bptree.child(ffIter, child);
      args.push_back(childIter->output * childIter->currWeight);
    }
    ffIter->output = ffIter->neuronTerminal->evaluate(args);
  }
  return bptree.begin()->output;

}


///
/// Creates back prop tree that will be used in running back propagation on 
/// this network.
/// @param postfix_stack Vector containing terminals which will be converted
/// to a more back propagation friendly format.
/// @param set Dataset
///
void BackPropAnnieTree::createTree(vector<TerminalSymbol*>& postfix_stack, Dataset* set){
  ExpressionTree extree;
  extree.convert_postfix(postfix_stack);
  
  bptree.clear();
  
  // first node in expression tree should be a neuron
  ExpressTreeIter exIter = extree.begin();
  Neuron topNeuron;
  topNeuron.neuronTerminal = exIter->el;
  bptree.insert(bptree.begin(), topNeuron);
  tree<Neuron>::iterator bpIter = bptree.begin();
  maxDepth = 0;

  addNeuron(exIter, bpIter, extree, bptree);

  // set adjusted depth for nodes based on maxDepth
  for(bpIter = bptree.begin(); bpIter != bptree.end(); bpIter++){
    bpIter->adjusted_depth = maxDepth - bptree.depth(bpIter);
  }
  
  // reset all leaves to be input layer (depth = 0)
  inputIndexes.clear();
  inputValue input;
  for(tree<Neuron>::leaf_iterator leafIter = bptree.begin_leaf(); leafIter != bptree.end_leaf();
    leafIter++){
    leafIter->adjusted_depth = 0;
    // go through leaves and extract symbols -- they will indicate which index to use
    // and whether Genotype ('G') or Covariate ('C')
    string inputString = leafIter->neuronTerminal->get_name();
    if(inputString[0] == 'G')
      input.isGenotype = true;
    else
      input.isGenotype = false;
    
    stringstream ss(inputString.substr(1, inputString.size()-1));
    ss >> input.index;
    input.index--; // index should be offset by one as it starts at zero
    inputIndexes.push_back(input);
  }
  
  // determine sizes of layers
  vector<int> layers(maxDepth+1,0);
  for(bpIter = bptree.begin(); bpIter != bptree.end(); bpIter++){
    layers[bpIter->adjusted_depth]++;
  }

  // construct multilayernetwork that matches the one created above
  // start with number of inputs
  if(net != NULL){
    delete net;
  }
  
  net = new MultiLayerNetwork(layers[0]);
  for(size_t currLayer = 1; currLayer < layers.size()-1; currLayer++){
    net->addLayer(layers[currLayer]);
  }
  // add single output node
  net->addLayer(1);
  
  // need to set the proper weights 
  // need indexes of neurons for connections
  vector<int> neuronIndexes(layers.size(), 0);
  vector<int> destIndexes(layers.size(), 0);
  
  // start at top and connect children to the root
  // for each child connect any children (and so on)
  // update indexes to match changes
  bpIter = bptree.begin();
  setMultiNetworkConnections(bpIter, bptree, neuronIndexes, destIndexes, net);
 
  constructTrainingSet(set);
}



///
/// Trains network using annie library
/// @param set Dataset
/// @param epochs Number of epochs to run
///
void BackPropAnnieTree::trainNetwork(Dataset* set, int epochs){ 
  // ready to train
try{
//  net->train(*trainSet, defaultControl["epochs"], defaultControl["learningRate"], defaultControl["momentum"]);'
  net->train(*trainSet, epochs, 0.3, 0.6);
}
catch(Exception &e){
  cout << e.what();
  cout << "Exit in trainNetwork" << endl;
  exit(1);
}
  
}

///
/// Constructs training set used by annie network
/// @set Dataset
///
void BackPropAnnieTree::constructTrainingSet(Dataset* set){
// cout << "construct training set " << trainSet << endl;
  if(trainSet != NULL){
    delete trainSet;
// cout << "deleting training set" << endl;
    trainSet=NULL;
  }
  // set size -- only one output in these networks
  trainSet = new TrainingSet(inputIndexes.size(), 1);
  Individual * currInd;
  
  for(unsigned int indindex=0; indindex < set->num_inds(); indindex++){
    currInd = set->get_ind(indindex);
    
    annie::Vector input;
    annie::Vector output;
    
    for(vector<inputValue>::iterator inputIter = inputIndexes.begin(); inputIter != inputIndexes.end();
      ++inputIter){
      
      if(inputIter->isGenotype)
        input.push_back(currInd->get_genotype(inputIter->index));
      else
        input.push_back(currInd->get_covariate(inputIter->index));
    }
    output.push_back(currInd->get_status());
    
    trainSet->addIOpair(input, output);
  }
}


///
/// Returns weights from current network
/// @return Weights in the back prop tree
///
vector<float> BackPropAnnieTree::getWeights(){
  // return weights from the MultiLayerNetwork
  // needs to be in same order as original weights
  vector<float> weights;
  
  tree<Neuron>::iterator bpIter = bptree.begin();
  ++bpIter; // skip output as it will not have any weights in the bptree representation
  
  Vector annieWeights;
  
  // retrieve weights
  for(;bpIter != bptree.end(); ++bpIter){
    Layer& layer = net->getLayer(bpIter->annie_layer);
    layer[bpIter->annie_neuron].getWeights(annieWeights);
    weights.push_back(annieWeights[bpIter->annie_link]);
    bpIter->currWeight = weights.back();
  }
  
  return weights;
}

///
/// Sets connections in MultiLayerNetwork recursively to
/// match original network
/// @param bpIter
/// @param bptree
/// @param sourceIndexes
/// @param destIndexes
/// @param net
///
void BackPropAnnieTree::setMultiNetworkConnections(tree<Neuron>::iterator& bpIter,
  tree<Neuron>& bptree, vector<int>& sourceIndexes, vector<int>& destIndexes,
  MultiLayerNetwork* net){

  tree<Neuron>::iterator childIter;
  int thisNeuronIndex = destIndexes[bpIter->adjusted_depth]++;
  
  int link=0;
  
  for(int child=0; child < int(bptree.number_of_children(bpIter)); child++){
    childIter = bptree.child(bpIter, child);
    
    // store info for retrieval after training back into bptree
    childIter->setAnnieLink(bpIter->adjusted_depth, thisNeuronIndex,
      link++);
    
    // connect child to current layer
    net->connect(childIter->adjusted_depth, sourceIndexes[childIter->adjusted_depth]++,
      bpIter->adjusted_depth, thisNeuronIndex, childIter->currWeight);
    
    // recursively connect if any children
    if(bptree.number_of_children(childIter) > 0){
      setMultiNetworkConnections(childIter, bptree, sourceIndexes, destIndexes, net);
    }
    
  }

}
