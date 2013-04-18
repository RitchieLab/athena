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
/// @param postFixStack Vector containing terminals which will be converted
/// into a more back propagation friendly format.
/// @param set Dataset
///
BackPropAnnieTree::BackPropAnnieTree(vector<TerminalSymbol*>& postFixStack, Dataset* set){
	initialize();
	createTree(postFixStack, set);
}


///
/// Initializes object
///
void BackPropAnnieTree::initialize(){
	beta = 0.3;
	maxDepth = 0;
	net = NULL;
	trainSet = NULL;
	biasValue = 1;
}


///
/// Recursively adds nodes to the backprop tree from the expression tree.
/// @param exIter iterator pointing to current neuron in the expression tree
/// @param bpIter iterator pointing to current neuron in backprop tree
/// @param exTree ExpressionTree 
/// @param bpTree Backpropagation tree
///
void BackPropAnnieTree::addNeuron(ExpressTreeIter& exIter, tree<Neuron>::iterator& bpIter,
	ExpressionTree& exTree, tree<Neuron>& bpTree){

	deque<float> args;

	for(int child=0; child < exTree.numberOfChildren(exIter); child++){
		Neuron newNode;
		
		ExpressTreeIter weightIter, childIter, exNodeIter;
		// each child will be a weight terminal symbol
		weightIter = exTree.child(exIter, child);
		bool any_children = false;
	 
		// now get the actual weight value and variable/neuron
		for(int wChild=0; wChild < exTree.numberOfChildren(weightIter); wChild++){
			childIter=exTree.child(weightIter, wChild);
			if(childIter->el->getTermType() == TerminalSymbol::Constant){
				newNode.constantWeightValue = childIter->el;
				newNode.currWeight = childIter->el->evaluate(args);
			}
			// in this case it is either a variable or another neuron
			else{
				newNode.neuronTerminal = childIter->el;
				exNodeIter = childIter;
				if(exTree.numberOfChildren(childIter) > 0)
					any_children=true;
			}
		}
		tree<Neuron>::iterator newBPnode = bpTree.append_child(bpIter, newNode);
		
		// check max depth here 
		if(bpTree.depth(newBPnode) > maxDepth){
			maxDepth = bpTree.depth(newBPnode);
		}
		
		// if any children of this neuron (it is a PADD) then recursively
		// call this function to continue building backprop tree
		if(any_children){
			addNeuron(exNodeIter, newBPnode, exTree, bpTree);
		}
	}
	
}


///
/// Calculates feed-forward result on individual specified
/// @param ind Individual whose values will be 
/// @return The output of network
///
float BackPropAnnieTree::feedForward(Individual* ind){
	
	ContinVariable::setInd(ind);
	GenotypeTerm::setInd(ind);
	
	tree<Neuron>::iterator childIter;
	deque<float> args;
	
	// need to begin at bottom of tree and then progress up
	// record output at each neuron
	for(tree<Neuron>::post_order_iterator ffIter=bpTree.begin_post();
		 ffIter != bpTree.end_post(); ffIter++){

		args.clear();
		
		// need to set output -- each neurons output is set by 
		// taking the outputs of neurons below -- adding together
		for(unsigned int child=0; child < bpTree.number_of_children(ffIter); child++){
			childIter = bpTree.child(ffIter, child);
			args.push_back(childIter->output * childIter->currWeight);
		}
		ffIter->output = ffIter->neuronTerminal->evaluate(args);
	}
	return bpTree.begin()->output;

}


///
/// Creates back prop tree that will be used in running back propagation on 
/// this network.
/// @param postFixStack Vector containing terminals which will be converted
/// to a more back propagation friendly format.
/// @param set Dataset
///
void BackPropAnnieTree::createTree(vector<TerminalSymbol*>& postFixStack, Dataset* set){
	ExpressionTree exTree;
	exTree.convertPostFix(postFixStack);
	
	bpTree.clear();
	
	// first node in expression tree should be a neuron
	ExpressTreeIter exIter = exTree.begin();
	Neuron topNeuron;
	topNeuron.neuronTerminal = exIter->el;
	bpTree.insert(bpTree.begin(), topNeuron);
	tree<Neuron>::iterator bpIter = bpTree.begin();
	maxDepth = 0;

	addNeuron(exIter, bpIter, exTree, bpTree);

	// set adjusted depth for nodes based on maxDepth
	for(bpIter = bpTree.begin(); bpIter != bpTree.end(); bpIter++){
		bpIter->adjustedDepth = maxDepth - bpTree.depth(bpIter);
	}
	
	// reset all leaves to be input layer (depth = 0)
	inputIndexes.clear();
	inputValue input;
	for(tree<Neuron>::leaf_iterator leafIter = bpTree.begin_leaf(); leafIter != bpTree.end_leaf();
		leafIter++){
		leafIter->adjustedDepth = 0;
		// go through leaves and extract symbols -- they will indicate which index to use
		// and whether Genotype ('G') or Covariate ('C')
		string inputString = leafIter->neuronTerminal->getName();
		if(inputString[0] == 'G'){
			input.inType = Genotype;
		}
		else if(inputString[0] == 'C'){
			input.inType = Contin;
		}
		else{
			input.inType = Bias;
			deque<float> args;
			biasValue = leafIter->neuronTerminal->evaluate(args);
		}

		stringstream ss(inputString.substr(1, inputString.size()-1));
		ss >> input.index;
		input.index--; // index should be offset by one as it starts at zero
		inputIndexes.push_back(input);
	}
	
	// determine sizes of layers
	vector<int> layers(maxDepth+1,0);
	for(bpIter = bpTree.begin(); bpIter != bpTree.end(); bpIter++){
		layers[bpIter->adjustedDepth]++;
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
	bpIter = bpTree.begin();
	setMultiNetworkConnections(bpIter, bpTree, neuronIndexes, destIndexes, net);
 
	constructTrainingSet(set);
}



///
/// Trains network using annie library
/// @param set Dataset
/// @param epochs Number of epochs to run
///
void BackPropAnnieTree::trainNetwork(Dataset* set, int epochs){ 
try{
	net->train(*trainSet, epochs, 0.3, 0.6);
}
catch(Exception &e){
	cerr << e.what();
	cerr << "Exit in trainNetwork" << endl;
	exit(1);
}
	
}

///
/// Constructs training set used by annie network
/// @set Dataset
///
void BackPropAnnieTree::constructTrainingSet(Dataset* set){
	if(trainSet != NULL){
		delete trainSet;
		trainSet=NULL;
	}
	// set size -- only one output in these networks
	trainSet = new TrainingSet(inputIndexes.size(), 1);
	Individual * currInd;
	
	for(unsigned int indIndex=0; indIndex < set->numInds(); indIndex++){
		currInd = set->getInd(indIndex);
		
		annie::Vector input;
		annie::Vector output;
		
		for(vector<inputValue>::iterator inputIter = inputIndexes.begin(); inputIter != inputIndexes.end();
			++inputIter){
			switch(inputIter->inType){
				case Genotype:
					input.push_back(currInd->getGenotype(inputIter->index));
					break;
				case Contin:
					input.push_back(currInd->getCovariate(inputIter->index));
					break;
				case Bias:
					input.push_back(biasValue);
					break;
			}  
				
		}
		output.push_back(currInd->getStatus());
		
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
	
	tree<Neuron>::iterator bpIter = bpTree.begin();
	// skip output as it will not have any weights in the bpTree representation
	++bpIter; 
	
	Vector annieWeights;
	
	// retrieve weights
	for(;bpIter != bpTree.end(); ++bpIter){
		Layer& layer = net->getLayer(bpIter->annieLayer);
		layer[bpIter->annieNeuron].getWeights(annieWeights);
		weights.push_back(annieWeights[bpIter->annieLink]);
		bpIter->currWeight = weights.back();
	}
	
	return weights;
}



///
/// Sets connections in MultiLayerNetwork recursively to
/// match original network
/// @param bpIter
/// @param bpTree
/// @param sourceIndexes
/// @param destIndexes
/// @param net
///
void BackPropAnnieTree::setMultiNetworkConnections(tree<Neuron>::iterator& bpIter,
	tree<Neuron>& bpTree, vector<int>& sourceIndexes, vector<int>& destIndexes,
	MultiLayerNetwork* net){

	tree<Neuron>::iterator childIter;
	int thisNeuronIndex = destIndexes[bpIter->adjustedDepth]++;
	
	int link=0;
	
	for(int child=0; child < int(bpTree.number_of_children(bpIter)); child++){
		childIter = bpTree.child(bpIter, child);
		
		// store info for retrieval after training back into bpTree
		childIter->setAnnieLink(bpIter->adjustedDepth, thisNeuronIndex,
			link++);
		
		// connect child to current layer
		net->connect(childIter->adjustedDepth, sourceIndexes[childIter->adjustedDepth]++,
			bpIter->adjustedDepth, thisNeuronIndex, childIter->currWeight);
		
		// recursively connect if any children
		if(bpTree.number_of_children(childIter) > 0){
			setMultiNetworkConnections(childIter, bpTree, sourceIndexes, destIndexes, net);
		}
	}
}

