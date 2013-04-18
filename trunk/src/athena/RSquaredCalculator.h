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
 * File:   RSquared.h
 * Author: dudeksm
 *
 */

#ifndef _RSQUAREDCALCULATOR_H
#define	_RSQUAREDCALCULATOR_H

#include "SolutionCalculator.h"


#include <iostream>

///
/// Calculates mean squared error
///
class RSquaredCalculator: public SolutionCalculator{
		
		public:
		
		RSquaredCalculator();
		
		/// resets calculator for new analysis set
		void reset();
		
		/// adds evaluation score to results
		void addIndScore(float score, float status);
		
		/// returns r-squared
		float getScore(){
				return 1-(squaredErrorTotal / ssTotal);
		}
		
		/// returns false so that the best is the smallest
		bool maxBest(){return true;}
		
		/// returns worst score
		float getWorst(){return -100;}
		
		/// Sets the SSTotal which is the square of the differences of the mean status and each ind status
		void setConstant(float constant){ssTotal = constant;}
		
private:
		
		float squaredErrorTotal, ssTotal;
		unsigned int totalIndsTested;
		
		
};

#endif	/* _MEANSQUAREDERRCALCULATOR_H */

