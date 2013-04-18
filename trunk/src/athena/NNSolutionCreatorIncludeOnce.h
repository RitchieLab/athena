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
// NNSolutionCreatorIncludeOnce.h

#ifndef _NNSOLUTIONCREATORINCLUDEONCE_H
#define	_NNSOLUTIONCREATORINCLUDEONCE_H

#include "NNSolutionCreator.h"

///
/// Inherits from the NNSolutionCreator class.  
/// Adds functionality allowing which assigns worst fitness to any network
/// that doesn't contain all the variables in the grammar.
///

class NNSolutionCreatorIncludeOnce: public NNSolutionCreator{

	public:
	
		NNSolutionCreatorIncludeOnce(){varsSet = false;}
	
		/// creates solution from vector of strings
		void establishSolution(vector<string>& symbols, Dataset* set);
		
		/// returns fitness score through evaluation of solution
		float evaluate(Dataset* set);
		
		/// passes vector of symbols that must be part of the neural network
		void requiredSymbols(vector<string>& symbols, Dataset* set);
		
		/// gets list of required symbols
		void restrict(vector<string>& vars){reqSymbols = vars;}
		
		
	private:
		
		vector<TerminalSymbol*> reqGenos, reqCovars;
		vector<string> reqSymbols;
		bool allVarsIncluded, varsSet;
		
};

#endif
