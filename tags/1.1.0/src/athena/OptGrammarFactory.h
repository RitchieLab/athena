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
 * Created on November 8, 2013, 3:00 PM
 */

#ifndef _OPTGRAMMARFACTORY_H
#define	_OPTGRAMMARFACTORY_H

#include "AthenaExcept.h"
#include <map>
#include <string>

class OptGrammar;
typedef OptGrammar* (createOptFunc)();

///
/// Returns calculator depending on string passed
///
class OptGrammarFactory{
		
public:
	const std::string& registerCalc(const std::string& key, createOptFunc* ptr);
	OptGrammar* create(const std::string& key);

	static OptGrammarFactory& getFactory(){static OptGrammarFactory f; return f;}
		
private:

		OptGrammarFactory(){}
		OptGrammarFactory(const OptGrammarFactory&);
		OptGrammarFactory& operator=(const OptGrammarFactory&);
		
		std::map<const std::string, createOptFunc*> creation_map;
};



#endif	/* _OPTGRAMMARFACTORY_H */

