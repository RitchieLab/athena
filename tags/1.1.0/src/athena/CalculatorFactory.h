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
 * File:   CalculatorFactory.h
 * Author: dudeksm
 *
 * Created on November 20, 2008, 4:00 PM
 */

#ifndef _CALCULATORFACTORY_H
#define	_CALCULATORFACTORY_H

#include "AthenaExcept.h"
#include <map>
#include <string>

class SolutionCalculator;
typedef SolutionCalculator* (createFunc)();

///
/// Returns calculator depending on string passed
///
class CalculatorFactory{
		
public:
	const std::string& registerCalc(const std::string& key, createFunc* ptr);
	SolutionCalculator* create(const std::string& key);

	static CalculatorFactory& getFactory(){static CalculatorFactory f; return f;}
		
private:

		CalculatorFactory(){}
		CalculatorFactory(const CalculatorFactory&);
		CalculatorFactory& operator=(const CalculatorFactory&);
		
		std::map<const std::string, createFunc*> creation_map;
};



#endif	/* _CALCULATORFACTORY_H */

