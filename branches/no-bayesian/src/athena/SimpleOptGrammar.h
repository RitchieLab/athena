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
/* 
 * File:   SimpleOptGrammar.h
 * Author: dudeksm
 *
 * Created on November 8, 2013, 3:46 PM
 */

#ifndef _SIMPLEOPTGRAMMAR_H
#define	_SIMPLEOPTGRAMMAR_H
#include "OptGrammar.h"
#include<deque>


class ConstantPointer{
	public: 
	
	ConstantPointer(){
		cons=NULL;
	}
	
	ConstantPointer(TerminalSymbol* c){
		cons=c;
		std::deque<float> blank;
		value = cons->evaluate(blank);
	}

	bool operator< (const ConstantPointer& other) const {
		return value < other.value;
  }
	TerminalSymbol* cons;
	float value;
};


///
/// Calculates balanced accuracy
///
class SimpleOptGrammar : public OptGrammarImp<SimpleOptGrammar>{
	 
public:
		
		SimpleOptGrammar();
		
		void terminalsFromConstant(float value, symbVector& optSymbols);
		
		void addNumberConstant(TerminalSymbol * sym);
		
private:
		static const string optGrammarName;
		std::set<ConstantPointer> constantSet;
		std::set<ConstantPointer>::iterator minIter;
		std::set<ConstantPointer>::reverse_iterator maxIter;
		optSymbol maxTerm, minTerm;
};


#endif	/* _SIMPLEOPTGRAMMAR_H */

