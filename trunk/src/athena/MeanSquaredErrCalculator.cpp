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
#include "MeanSquaredErrCalculator.h"

using namespace std;

MeanSquaredErrCalculator::MeanSquaredErrCalculator(){
		reset();
		name = "Mean Squared Error";
}


void MeanSquaredErrCalculator::reset(){
	 totalIndsTested = 0;
	 squaredErrorTotal = 0.0;
	 ssTotal = 0.0;
	 statTotal.clear();
}


///
/// Adds score to running total within object
/// @param score
///
void MeanSquaredErrCalculator::addIndScore(float score, float stat){
		float difference = score - stat;
		squaredErrorTotal += difference * difference;
		totalIndsTested++;
		ssTotal += stat;
		statTotal.push_back(stat);
}


///
/// Calculates and returns ssTotal for use in calculating
/// R-squared value
/// @return ssTotal
///
float MeanSquaredErrCalculator::getConstant(){
	float mean = ssTotal / totalIndsTested;
	float diff=0.0;
	for(vector<float>::iterator iter=statTotal.begin(); iter != statTotal.end();
		++iter){
		diff = diff + (*iter-mean) * (*iter-mean);
	}  
	ssTotal = diff;
	return ssTotal;
}
