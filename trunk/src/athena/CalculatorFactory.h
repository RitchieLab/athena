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

#include "SolutionCalculator.h"
#include "AthenaExcept.h"
#include <map>
#include <string>

///
/// Returns calculator depending on string passed
///
class CalculatorFactory{
		
public:
				
		static SolutionCalculator* createCalculator(std::string calcName);
		
private:
		
		static void setCalcMap();
		
		enum CalcType{
				NoCalcType,
				BalanceCalcType,
				MeanSquaredErrType,
				RSquaredType
		};
		
		static std::map<std::string, CalcType> CalcMap;

};



#endif	/* _CALCULATORFACTORY_H */

