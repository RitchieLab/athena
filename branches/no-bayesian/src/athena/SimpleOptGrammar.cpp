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
#include "SimpleOptGrammar.h"
#include<sstream>
#include<iostream>

string const SimpleOptGrammar::optGrammarName = SimpleOptGrammar::registerCalc("SIMPLE");


SimpleOptGrammar::SimpleOptGrammar(){
	maxTerm.noNT = false;
	minTerm.noNT = false;
}


void SimpleOptGrammar::addNumberConstant(TerminalSymbol * sym){
	constantSet.insert(sym);
	maxIter = constantSet.rbegin();
	minIter = constantSet.begin();
	maxTerm.symbol = maxIter->cons->getName();
	minTerm.symbol = minIter->cons->getName();
}


void SimpleOptGrammar::terminalsFromConstant(float value, 
	symbVector& optSymbols){
	
	optSymbols.clear();
	
	if(value >= maxIter->value){
		optSymbols.push_back(maxTerm);
		return;
	}
	else if(value <= minIter->value){
		optSymbols.push_back(minTerm);
		return;
	}
	
	optSymbol constantTerm;
	// need to confirm this
	constantTerm.noNT = false;
	
	// get the closest constant to the value and place that symbol
	ConstantPointer temp;
	temp.value = value;
	std::set<ConstantPointer>::iterator upperIter = constantSet.upper_bound(temp);
	std::set<ConstantPointer>::reverse_iterator lowerIter = set<ConstantPointer>::reverse_iterator(upperIter);

	if((value-lowerIter->value) < (upperIter->value-value)){
		constantTerm.symbol = lowerIter->cons->getName();
	}
	else{
		constantTerm.symbol = upperIter->cons->getName();
	}
	optSymbols.push_back(constantTerm);
}	
