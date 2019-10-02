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
#include "K2Calculator.h"
#include<sstream>

string const K2Calculator::calcMatchName = K2Calculator::registerCalc("K2");


K2Calculator::K2Calculator(){
		reset();
		name = "K2";
		outputNames.push_back("not-improved");
// 		outputNames.push_back("bal-acc");
}


void K2Calculator::reset(){
	k2Score = 0.0;
}


///
/// Returns vector containing any additional values for output
///
std::vector<std::string> K2Calculator::getAdditionalFinalOutput(){
	return outputValues;
}



///
/// K2 currently doesn't have additional output to report
/// @param results stat::TestResult
///
void K2Calculator::evaluateAdditionalOutput(std::vector<stat::TestResult>& results){
	outputValues.clear();
// 	float auc = stat::AUCCalc::calculateAUC(results);
// 	std::stringstream ss;
// 	ss << auc;
// 	outputValues.push_back(ss.str());
}

///
/// Adds score to running total within object
/// @param score
///
float K2Calculator::addIndScore(float score, float stat){
	 k2Score += score;
	 return score;
}

