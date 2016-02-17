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
 * File:   SolutionCalculator.h
 * Author: dudeksm
 *
 * Created on November 20, 2008, 2:38 PM
 */

#ifndef _SOLUTIONCALCULATOR_H
#define	_SOLUTIONCALCULATOR_H

#include<string>
#include<vector>
#include <AUCcalc.h>
#include "CalculatorFactory.h"
#include <Dataset.h>


///
/// Base class for calculation of solution scores.
/// For example, calculates Balanced accuracy for binary outcomes and
/// mean squared error for continuous outcomes.
///

class SolutionCalculator{

public:

		SolutionCalculator(){name = "Calculator";}

		virtual ~SolutionCalculator(){}

		virtual void addIndScore(float score, float status)=0;

		virtual float getScore()=0;

		virtual void reset()=0;

		virtual bool maxBest()=0;

		virtual bool logMaxBest()=0;

		virtual float getWorst()=0;

		virtual float getConstant(){return 0.0;}

		virtual bool requiresCaseControl(){return false;}

		virtual SolutionCalculator* clone()=0;

		/// Used when a sub class needs a constant value for calculations as in RSquared
		virtual void setConstant(data_manage::Dataset* ds)=0;

		virtual std::vector<std::string> getAdditionalOutputNames(){
			return outputNames;
		}

		virtual std::vector<std::string> getAdditionalFinalOutput(){
			return outputValues;
		}

		virtual void evaluateAdditionalOutput(std::vector<stat::TestResult>& results){}

		std::string getName(){return name;}

protected:
		std::string name;
		std::vector<std::string> outputNames, outputValues;

};


template <class T>
class SolutionCalcImp : public SolutionCalculator {
public:
	static SolutionCalculator* create(){return new T();}

protected:
	static const std::string& registerCalc(const std::string& key_in);
};

template<typename T>
const std::string& SolutionCalcImp<T>::registerCalc(const std::string& keyIn){
	return CalculatorFactory::getFactory().registerCalc(keyIn, &T::create);
}

#endif	/* _SOLUTIONCALCULATOR_H */

