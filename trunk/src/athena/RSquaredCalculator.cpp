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
#include "RSquaredCalculator.h"

RSquaredCalculator::RSquaredCalculator(){
		reset();
		ssTotal = 1;
		name = "R-Squared";
}

void RSquaredCalculator::reset(){
	 totalIndsTested = 0;
	 squaredErrorTotal = 0.0;
}

///
/// Adds score to running total within object
/// @param score
///
void RSquaredCalculator::addIndScore(float score, float stat){
		float difference = score - stat;
		squaredErrorTotal += difference * difference;
		totalIndsTested++;
}
