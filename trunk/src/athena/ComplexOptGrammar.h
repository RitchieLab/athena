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
 * File:   ComplexOptGrammar.h
 * Author: dudeksm
 *
 * Created on November 8, 2013, 3:46 PM
 */

#ifndef _COMPLEXOPTGRAMMAR_H
#define	_COMPLEXOPTGRAMMAR_H
#include "OptGrammar.h"
#include "Structs.h"

///
/// Calculates balanced accuracy
///
class ComplexOptGrammar : public OptGrammarImp<ComplexOptGrammar>{
	 
public:
		
		ComplexOptGrammar();
		
		void terminalsFromConstant(float value, symbVector& optSymbols);
		
private:
		void getStringFromNum(float value, symbVector& optSymbols);

		static const string optGrammarName;
		optSymbol leftParenSymb, rightParenSymb, concatSymb, periodSymb;
};


#endif	/* _BALACCCALCULATOR_H */

