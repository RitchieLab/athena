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
#include "K2AICCalculator.h"
#include<sstream>

string const K2AICCalculator::calcMatchName = K2AICCalculator::registerCalc("K2-AIC");


K2AICCalculator::K2AICCalculator(){
		reset();
		name = "K2-AIC";
		kTerm = 1.0;
		outputNames.push_back("not-improved");
		outputNames.push_back("bal-acc");
}



void K2AICCalculator::reset(){
	k2Score = 0.0;
}


///
/// Returns vector containing any additional values for output
///
std::vector<std::string> K2AICCalculator::getAdditionalFinalOutput(){
	return outputValues;
}



///
/// K2 currently doesn't have additional output to report
/// @param results stat::TestResult
///
void K2AICCalculator::evaluateAdditionalOutput(std::vector<stat::TestResult>& results){
	outputValues.clear();
}

///
/// Adds score to running total within object
/// @param score
/// @param stat is number of parameters to use in penalty term
///
void K2AICCalculator::addIndScore(float score, float stat){
	 k2Score += score;
	 // penalty term
	 k2Score -= kTerm * stat;
}

