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
 * File:   TerminalSymbCreator.h
 * Author: dudeksm
 *
 * Created on November 19, 2008, 4:00 PM
 */

#ifndef _TERMINALSYMBCREATOR_H
#define	_TERMINALSYMBCREATOR_H

#include "TerminalSymbol.h"
#include <map>
#include <vector>
#include "AthenaExcept.h"
#include <Individual.h>
#include "Structs.h"

using namespace data_manage;

class TerminalSymbCreator{
		
	public:
		TerminalSymbCreator();
		void createTerminals(int numGenotypes, int numCovariates);
		TerminalSymbol * createConstant(const std::string& symbol);
		
		void terminalsFromConstant(float value, symbVector& optSymbols);
		
		void addGenotypeVariables(int numVariables);
		void addContinVariables(int numVariables);
		
		TerminalSymbol* getTerm(std::string& symbol);
		
		inline TerminalSymbol* rightParen(){return rParen;}
		inline TerminalSymbol* leftParen(){return lParen;}
		inline TerminalSymbol* comma(){return commaPtr;}
		inline TerminalSymbol* concaten(){return concat;}
		
		void setInd(Individual* ind);
		
	private:
	 
		void getStringFromNum(float value, symbVector& optSymbols);
		std::map<std::string, TerminalSymbol*> terminalMap;
		TerminalSymbol* rParen, *lParen, *commaPtr, *concat;
		optSymbol leftParenSymb, rightParenSymb, concatSymb, periodSymb;

};



#endif	/* _TERMINALSYMBCREATOR_H */

