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
 * File:   OptGrammar.h
 * Author: dudeksm
 *
 * Created on November 20, 2008, 2:38 PM
 */

#ifndef _OPTGRAMMAR_H
#define	_OPTGRAMMAR_H

#include<string>
#include<vector>
#include "OptGrammarFactory.h"
#include "Structs.h"

#include "Terminals.h"

///
/// Base class for creating grammar from optimizing weights
///

class OptGrammar{
		
public:
		
		OptGrammar(){}
		
		virtual ~OptGrammar(){}

		virtual void terminalsFromConstant(float value, symbVector& optSymbols)=0;
		
		virtual void addNumberConstant(TerminalSymbol* sym){}
		
protected:

		
};


template <class T>
class OptGrammarImp : public OptGrammar {
public:
	static OptGrammar* create(){return new T();}

protected:
	static const std::string& registerCalc(const std::string& key_in);
};

template<typename T>
const std::string& OptGrammarImp<T>::registerCalc(const std::string& keyIn){
	return OptGrammarFactory::getFactory().registerCalc(keyIn, &T::create);
}

#endif	/* _OPTGRAMMAR_H */

