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
 * File:   MeanSquaredErrCalculator.h
 * Author: dudeksm
 *
 * Created on November 20, 2008, 3:45 PM
 */

#ifndef _MEANSQUAREDERRCALCULATOR_H
#define	_MEANSQUAREDERRCALCULATOR_H

#include "SolutionCalculator.h"


#include <iostream>
#include <vector>

///
/// Calculates mean squared error
///
class MeanSquaredErrCalculator: public SolutionCalcImp<MeanSquaredErrCalculator>{
		
		public:
		
		MeanSquaredErrCalculator();
		
		/// resets calculator for new analysis set
		void reset();
		
		/// adds evaluation score to results
		void addIndScore(float score, float status);
		
		/// returns mean squared error
		float getScore(){
				if(totalIndsTested == 0)
					return getWorst();
				else
					return squaredErrorTotal / totalIndsTested;
		}
		
		/// returns false so that the best is the smallest
		bool maxBest(){return false;}
		
		bool logMaxBest(){return true;}
		
		/// returns worst score
		float getWorst(){return 1000000;}
		/// sets the sstotal for use in calculating rsquared
		void setConstant(data_manage::Dataset* ds){ssTotal = ds->getSSTotal();}
		/// returns the sstotal for use in calculating rsquared
		float getConstant();

		/// returns rsquared for display to user
		float outputScore(float sc){return 1-(sc/ssTotal);}

		
private:

		static const string calcMatchName;
		
		float squaredErrorTotal, ssTotal;
		std::vector<float> statTotal;
		unsigned int totalIndsTested;
		
		
};




#endif	/* _MEANSQUAREDERRCALCULATOR_H */

