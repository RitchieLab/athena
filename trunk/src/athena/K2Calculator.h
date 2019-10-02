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
 * File:   K2Calculator.h
 * Author: dudeksm
 *
 * Created on May, 7 2014, 11:10 AM
 */

#ifndef _K2CALCULATOR_H
#define	_K2CALCULATOR_H

#include "SolutionCalculator.h"

///
/// Calculates balanced accuracy
///
class K2Calculator : public SolutionCalcImp<K2Calculator>{

public:

		K2Calculator();

		/// resets calculator for new analysis set
		void reset();

		/// adds evaluation score to results
		float addIndScore(float score, float status);

		/// returns balanced accuracy
		float getScore(){
			return -k2Score;
		}

	 /// scores are converted to positive values so max is best
	 bool maxBest(){return true;}

	 K2Calculator* clone(){return new K2Calculator(*this);}

	 bool logMaxBest(){return true;}

	 	virtual std::vector<std::string> getAdditionalOutputNames(){
			return outputNames;
		}

		virtual std::vector<std::string> getAdditionalFinalOutput();

	 	virtual void evaluateAdditionalOutput(std::vector<stat::TestResult>& results);

	 	virtual bool requiresCaseControl(){return false;}

	 	void setConstant(data_manage::Dataset* ds){}

	 /// returns worst score
	 float getWorst(){return 0.0;}
private:
		static const string calcMatchName;
		double k2Score;
};


#endif	/* _K2CALCULATOR_H */

