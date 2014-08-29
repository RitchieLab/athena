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
#include "K2BICCalculator.h"
#include<sstream>

string const K2BICCalculator::calcMatchName = K2BICCalculator::registerCalc("K2-BIC");


K2BICCalculator::K2BICCalculator(){
		reset();
		name = "K2-BIC";
		kTerm = 1.0;
		outputNames.push_back("not-improved");
		outputNames.push_back("bal-acc");
}



void K2BICCalculator::reset(){
	k2Score = 0.0;
}


///
/// Returns vector containing any additional values for output
///
std::vector<std::string> K2BICCalculator::getAdditionalFinalOutput(){
	return outputValues;
}



///
/// K2 currently doesn't have additional output to report
/// @param results stat::TestResult
///
void K2BICCalculator::evaluateAdditionalOutput(std::vector<stat::TestResult>& results){
	outputValues.clear();
}

///
/// Adds score to running total within object
/// @param score
/// @param stat is number of parameters to use in penalty term
///
void K2BICCalculator::addIndScore(float score, float stat){
	 k2Score += score;
	 // penalty term
	 k2Score -= kTerm * stat;
}

