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
//NNSolutionCreatorIncludeAll.cpp

#include <iostream>
using namespace std;

#include "NNSolutionCreatorIncludeAll.h"

///
/// Establishes solution.  If not all variables are in set, doesn't bother
/// evaluating model at all.  Sets flag so that worst fitness will be returned.
///
void NNSolutionCreatorIncludeAll::establishSolution(vector<string>& symbols, Dataset* set){
	
	NNSolutionCreator::establishSolution(symbols, set);
		
	if(!varsSet){
		requiredSymbols(reqSymbols, set);
		varsSet = true;
	}

	// check genos
	allVarsIncluded = true;
	vector<TerminalSymbol*>::iterator iter;
	for(iter=reqGenos.begin(); iter != reqGenos.end(); *iter++){
		if(genos.find(*iter) == genos.end()){
			allVarsIncluded = false;
			return;
		}
	}
	
	// check covars
	for(iter=reqCovars.begin(); iter != reqCovars.end(); *iter++){
		if(covars.find(*iter) == covars.end()){
			allVarsIncluded = false;
			return;
		}
	}
}

///
/// returns fitness score through evaluation of solution
///
float NNSolutionCreatorIncludeAll::evaluate(Dataset* set){

	if(allVarsIncluded)
		return NNSolutionCreator::evaluate(set);
	else
		return calculator->getWorst();
}


///
/// Takes list of symbols and Dataset to establish list of terminals that must 
/// be part of the neural network for it to receive 
///
void NNSolutionCreatorIncludeAll::requiredSymbols(vector<string>& symbols, Dataset* set){
	
	if(!terminalsSet)
		NNSolutionCreator::setVariables(set);
		
	// use symbols to get pointers that can be used to check for inclusion of all variables
	reqCovars.clear();
	reqGenos.clear();
	vector<string>::iterator iter;
	TerminalSymbol* var;
	for(iter=symbols.begin(); iter != symbols.end(); iter++){
		var = termHolder.getTerm(*iter);
		if(var->getTermType() == TerminalSymbol::Genotype)
			reqGenos.push_back(var);
		else if(var->getTermType() == TerminalSymbol::Covariate)
			reqCovars.push_back(var);
	}
}
