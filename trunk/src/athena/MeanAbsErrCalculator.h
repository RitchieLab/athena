/*
Copyright Marylyn Ritchie 2014

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
 * File:   MeanAbsErrCalculator.h
 * Author: dudeksm
 *
 * Created on June 20, 2014, 3:45 PM
 */

#ifndef _MEANABSERRCALCULATOR_H
#define	_MEANABSERRCALCULATOR_H

#include "SolutionCalculator.h"


#include <iostream>
#include <vector>

///
/// Calculates mean squared error
///
class MeanAbsErrCalculator: public SolutionCalcImp<MeanAbsErrCalculator>{
		
		public:
		
		MeanAbsErrCalculator();
		
		/// resets calculator for new analysis set
		void reset();
		
		/// adds evaluation score to results
		void addIndScore(float score, float status);
		
		/// returns mean squared error
		float getScore(){
				if(totalIndsTested == 0)
					return getWorst();
				else
					return  errTotal / obsTotal;
		}
		
		/// returns false so that the best is the smallest
		bool maxBest(){return false;}
		
		bool logMaxBest(){return true;}
		
		/// returns worst score
		float getWorst(){return 1e6;}
		/// sets the sstotal for use in calculating rsquared
		void setConstant(data_manage::Dataset* ds){}
		/// returns the sstotal for use in calculating rsquared
		float getConstant(){ return 0.0;}

		/// returns 1-score as report value
		float outputScore(float sc){return 1-sc;}

		
private:

		static const string calcMatchName;	
		float errTotal, obsTotal;
		unsigned int totalIndsTested;
};




#endif	/* _MEANABSERRCALCULATOR_H_ */

