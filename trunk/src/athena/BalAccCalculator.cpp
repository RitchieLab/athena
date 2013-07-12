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
#include "BalAccCalculator.h"
#include<sstream>

string const BalAccCalculator::calcMatchName = BalAccCalculator::registerCalc("BALANCEDACC");


BalAccCalculator::BalAccCalculator(){
		reset();
		name = "AUC";
		outputNames.push_back("AUC");
}



void BalAccCalculator::reset(){
		caseRight=0;
		caseWrong=0;
		controlRight=0;
		controlWrong=0;
}


///
/// Returns vector containing any additional values for output
///
std::vector<std::string> BalAccCalculator::getAdditionalFinalOutput(){
	return outputValues;
}



///
/// Calculates AUC and stores formatted output
/// @param results stat::TestResult
///
void BalAccCalculator::evaluateAdditionalOutput(std::vector<stat::TestResult>& results){
	outputValues.clear();
	float auc = stat::AUCCalc::calculateAUC(results);
	std::stringstream ss;
	ss << auc;
	outputValues.push_back(ss.str());
}

///
/// Adds score to running total within object
/// @param score
///
void BalAccCalculator::addIndScore(float score, float stat){
	 
		unsigned int result = score > 0.5?1:0;
		unsigned int status = (unsigned int)stat;
		
		if(result != status){
				if(status)
						caseWrong++;
				else
						controlWrong++;
		}
		else if(status)
				caseRight++;
		else
				controlRight++;
}

