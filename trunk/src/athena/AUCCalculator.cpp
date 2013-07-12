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
#include "AUCCalculator.h"

AUCCalculator::AUCCalculator(){
		reset();
		name = "AUC";
		outputNames.push_back("AUC");
}

void AUCCalculator::reset(){
	results.clear();
	auc = 0.0;
}

///
/// Calculates AUC and returns score
/// @return AUC
///
float AUCCalculator::getScore(){
	float auc = stat::AUCCalc::calculateAUC(results);
	return auc;
}

///
/// Adds score to running total within object
/// @param score for the individual evaluation
/// @param stat Status for the individual evaluation
///
void AUCCalculator::addIndScore(float score, float stat){ 
	stat::TestResult tempResult;
	tempResult.score = score;
	tempResult.status = stat;
	results.push_back(tempResult);
}

