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
 * File:   BalAccCalculator.h
 * Author: dudeksm
 *
 * Created on November 20, 2008, 2:46 PM
 */

#ifndef _BALACCCALCULATOR_H
#define	_BALACCCALCULATOR_H

#include "SolutionCalculator.h"

///
/// Calculates balanced accuracy
///
class BalAccCalculator : public SolutionCalcImp<BalAccCalculator>{

public:

		BalAccCalculator();

		/// resets calculator for new analysis set
		void reset();

		/// adds evaluation score to results
		float addIndScore(float score, float status);

		/// returns balanced accuracy
		float getScore(){
			return .5 * (float(caseRight) / (caseRight + caseWrong) +
			float(controlRight) /(controlRight + controlWrong));
		}
	 bool maxBest(){return true;}

	 bool logMaxBest(){return true;}

	 BalAccCalculator* clone(){return new BalAccCalculator(*this);}

	 	virtual std::vector<std::string> getAdditionalOutputNames(){
			return outputNames;
		}

		virtual std::vector<std::string> getAdditionalFinalOutput();

	 	virtual void evaluateAdditionalOutput(std::vector<stat::TestResult>& results);

	 	virtual bool requiresCaseControl(){return true;}

	 	void setConstant(data_manage::Dataset* ds){}

	 /// returns worst score
	 float getWorst(){return 0.0;}
private:
		static const string calcMatchName;
		float balancedAcc;
		unsigned int caseRight, caseWrong, controlRight, controlWrong;
};


#endif	/* _BALACCCALCULATOR_H */

