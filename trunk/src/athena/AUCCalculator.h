/*
Copyright Marylyn Ritchie 2013

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
 * File:   AUCCalculator.h
 * Author: dudeksm
 *
 */

#ifndef _AUCCALCULATOR_H
#define	_AUCCALCULATOR_H

#include "SolutionCalculator.h"

///
/// Calculates balanced accuracy
///
class AUCCalculator: public SolutionCalcImp<AUCCalculator>{
	 
public:
		
		AUCCalculator();
		
		/// resets calculator for new analysis set
		void reset();
		
		/// adds evaluation score to results
		void addIndScore(float score, float status);
		
		/// returns area under the curve
		float getScore();
		
	 	bool maxBest(){return true;}
		
	 /// returns worst score
	 float getWorst(){return 0.0;}
	 
	 void evaluateAdditionalOutput(std::vector<stat::TestResult>& results);
	 
	 virtual bool requiresCaseControl(){return true;}
	 
private:
		float auc;
		static const string calcMatchName;
		std::vector<stat::TestResult> results;
};


#endif	/* _BALACCCALCULATOR_H */

