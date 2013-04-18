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
// BackPropAnnie.h
#ifndef __BACKPROPANNIE_H__
#define __BACKPROPANNIE_H__

#include "Terminals.h"
#include <Dataset.h>
//#include "BackPropTree.h"
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
		int runBackProp(vector<TerminalSymbol*> & postFixStack, Dataset * set);

		/// Returns the list of weights in network
		vector<float> getWeights(){return bpAnnieTree.getWeights();}
		
		/// Returns optimized score
		float getOptimizedScore(){return optimizedScore;}

	private:
		
		void initialize();
		
		BackPropAnnieTree bpAnnieTree;

		double calculateMSE(Dataset* set);
		float optimizedScore;
};

#endif
