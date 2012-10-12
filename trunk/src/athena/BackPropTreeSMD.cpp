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
#include "BackPropTree.h"
#include <deque>

BackPropTree::BackPropTree(){
  initialize();
}


///
/// Alternative constructor.  Creates a back prop tree from the 
/// postfix stack.
/// @param postfix_stack Vector containing terminals which will be converted
/// into a more back propagation friendly format.
///
BackPropTree::BackPropTree(vector<TerminalSymbol*>& postfix_stack){
  initialize();
  createTree(postfix_stack);
}


///
/// Initializes object
///
void BackPropTree::initialize(){
  beta = 0.5;
}


///
/// Recursively adds nodes to the backprop tree from the expression tree.
/// @param exIter iterator pointing to current neuron in the expression tree
/// @param bpIter iterator pointing to current neuron in backprop tree
/// @param extree ExpressionTree 
/// @param bptree Backpropagation tree
///
void BackPropTree::addNeuron(ExpressTreeIter& exIter, tree<Neuron>::iterator& bpIter,
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
float BackPropTree::feedForward(Individual* ind){
  
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
/// Propagates errors down through tree and adjusts weights
/// @param ind Individual to use for back propagation
/// @param Return the output of the tree
///
float BackPropTree::backPropagate(Individual* ind){
  
  float output = feedForward(ind);
  tree<Neuron>::iterator bpIter = bptree.begin();
  
  // find delta for output layer
  bpIter->delta = output * (1-output)*(ind->status()-output);
  
  tree<Neuron>::iterator parentIter;
  
  // find delta for layers below
  for(++bpIter; bpIter != bptree.end(); bpIter++){
    float sum = 0.0;
    parentIter = bptree.parent(bpIter);

    //Only calculate the delta if node has any children 
    if(bptree.number_of_children(bpIter) > 0){
    // sum is the value of weights going up tree * delta of the node above
      sum = parentIter->delta * bpIter->currWeight;    
      bpIter->delta = bpIter->output * (1-bpIter->output) * sum;   
    }
  }

  // stop the loop when reach top node as it doesn't have any weights 
  // to update -- that is why the postIter is being checked against the 
  // standard begin function
  for(tree<Neuron>::post_order_iterator postIter=bptree.begin_post();
     postIter != bptree.begin(); postIter++){
     // the weights need to be evaluated here based on the deltas 
     // and previous weights
     parentIter = bptree.parent(postIter);
     postIter->prevWeight = beta * parentIter->delta * postIter->output;
     postIter->currWeight += postIter->prevWeight;
  }
  return output;
}

///
/// Returns weights from current network
/// @return Weights in the back prop tree
///
vector<float> BackPropTree::getWeights(){
  
  vector<float> weights;
  tree<Neuron>::iterator iter=bptree.begin();

  // top node has no weights as weights point up to higher nodes
  for(++iter; iter != bptree.end();
    ++iter){
    weights.push_back(iter->currWeight);
  }
  
  return weights;
}


///
/// Creates back prop tree that will be used in running back propagation on 
/// this network.
/// @param postfix_stack Vector containing terminals which will be converted
/// to a more back propagation friendly format.
///
void BackPropTree::createTree(vector<TerminalSymbol*>& postfix_stack){
  ExpressionTree extree;
  extree.convert_postfix(postfix_stack);
  
  bptree.clear();
  
  // first node in expression tree should be a neuron
  ExpressTreeIter exIter = extree.begin();
  Neuron topNeuron;
  topNeuron.neuronTerminal = exIter->el;
  bptree.insert(bptree.begin(), topNeuron);
  tree<Neuron>::iterator bpIter = bptree.begin();
  
  addNeuron(exIter, bpIter, extree, bptree);

}

