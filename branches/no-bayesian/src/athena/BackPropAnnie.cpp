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
/// @param postFixStack vector of TerminalSymbols representing neural network
/// @param set Dataset used for running backpropagation
/// @returns number of epochs
///
int BackPropAnnie::runBackProp(vector<TerminalSymbol*> & postFixStack, Dataset* set){

	// create the back propagation tree that will be used to adjust the weights 
	// in the original neural network
	bpAnnieTree.createTree(postFixStack, set);
	
	float mse, lastMSE=1000.0;

	int numEpochs=0;
	
	for(int i=0; i<100; i++){
		bpAnnieTree.trainNetwork(set, 1);
		numEpochs += 1;
		bpAnnieTree.getWeights();
		mse = calculateMSE(set);
		// original difference was 0.0001
		// stop when Mean-squared error is not changing or the mean-squared error is getting bigger
		if(fabs(mse-lastMSE)<threshold || mse > lastMSE){
			break;
		}
		else{
			lastMSE = mse;
		}
	}  
	
	optimizedScore = mse;
	return numEpochs;
}



///
/// Calculates mean-squared error on current backprop tree using 
/// current weights.
/// @param set Dataset used for running back prop
/// @return mse
///
double BackPropAnnie::calculateMSE(Dataset* set){
	double total=0.0, output;
	for(unsigned int indIndex=0; indIndex < set->numInds(); indIndex++){
		output =  bpAnnieTree.feedForward((*set)[indIndex]);
		total += (output - (*set)[indIndex]->status()) * (output - (*set)[indIndex]->status());
	}
	return  total/set->numInds();
	
}


