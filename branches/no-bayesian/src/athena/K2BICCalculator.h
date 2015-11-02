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
 * File:   K2BICCalculator.h
 * Author: dudeksm
 *
 * Created on May, 15 2014, 4:10 PM
 */

#ifndef _K2BICCalculator_H
#define	_K2BICCalculator_H

#include "SolutionCalculator.h"
#include <math.h>

///
/// Calculates balanced accuracy
///
class K2BICCalculator : public SolutionCalcImp<K2BICCalculator>{
	 
public:
		
		K2BICCalculator();
		
		/// resets calculator for new analysis set
		void reset();
		
		/// adds evaluation score to results
		void addIndScore(float score, float status);
		
		/// returns K2 with BIC penalty
		float getScore(){
			return -k2Score;
		}
		
	 /// K2 scores are negative so better are more negative? or less negative?
	 bool maxBest(){return true;}
	 
	 bool logMaxBest(){return true;}
	 
	 	virtual std::vector<std::string> getAdditionalOutputNames(){
			return outputNames;
		}
		
		virtual std::vector<std::string> getAdditionalFinalOutput();

	 	virtual void evaluateAdditionalOutput(std::vector<stat::TestResult>& results);
	 	
	 	virtual bool requiresCaseControl(){return false;}
	 	
	 	void setConstant(data_manage::Dataset* ds){kTerm=log(double(ds->numInds()))/2;}
	 
	 /// returns worst score
	 float getWorst(){return 0.0;}
private:
		static const string calcMatchName;
		double k2Score, kTerm;
};


#endif	/* _K2BICCalculator_H */

