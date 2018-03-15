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
#include "OptGrammar.h"

using namespace data_manage;


class TerminalSymbCreator{

	public:
		TerminalSymbCreator();
		~TerminalSymbCreator();
		void createTerminals(int numGenotypes, int numCovariates);
		TerminalSymbol * createConstant(const std::string& symbol);
		void addConstant(const string& symbol);

		void terminalsFromConstant(float value, symbVector& optSymbols);

		void addGenotypeVariables(int numVariables);
		void addContinVariables(int numVariables);

		TerminalSymbol* getTerm(std::string& symbol);

		inline TerminalSymbol* rightParen(){return rParen;}
		inline TerminalSymbol* leftParen(){return lParen;}
		inline TerminalSymbol* comma(){return commaPtr;}
		inline TerminalSymbol* concaten(){return concat;}
		inline TerminalSymbol* connector(){return connect;}
		inline TerminalSymbol* phenotype(){return phenoPtr;}

		std::string getGenoName(int varIndex);
		std::string getContinName(int varIndex);

		void setInd(Individual* ind);

// 		TerminalSymbol* getClosestConstant(float value);

		void setGrammerOptimization(string optName);

	private:

		std::map<std::string, TerminalSymbol*> terminalMap;
		TerminalSymbol* rParen, *lParen, *commaPtr, *concat, *connect, *phenoPtr;
		OptGrammar* grammarOptimizer;


};

#endif	/* _TERMINALSYMBCREATOR_H */

