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
//BackProp.h

#ifndef __BACKPROP_H__
#define __BACKPROP_H__

#include "Terminals.h"
#include <Dataset.h>

class BackProp{

	public:
	
		BackProp();
	
		virtual ~BackProp(){}
	
		/// Receives a postfix stack representation of the neural network
		virtual int runBackProp(vector<TerminalSymbol*> & postFixStack, Dataset * set)=0;
	
		/// Sets learning rate for the network
		void setLearningRate(double newBeta){learnRate = newBeta;}
	
		/// Returns learning rate
		double getLearningRate(){return learnRate;}
		
		/// Sets the number of iterations for running back propagation
		void setNumEpochs(int nEpochs){numEpochs = nEpochs;}
		
		/// Returns the number of iterations
		int getNumEpochs(){return numEpochs;}
		
		/// Sets threshold for stopping back propagation early
		void setBPThreshold(double newThresh){threshold = newThresh;}
		
		/// Returns the threshold
		double getBPThreshold(){return threshold;}
		
		/// Returns the list of weights in network
		virtual vector<float> getWeights()=0;    

	protected:
		void initialize();
		virtual double calculateMSE(Dataset* set)=0;

		double threshold, learnRate;
		int numEpochs, stochasticSize;
};

#endif
