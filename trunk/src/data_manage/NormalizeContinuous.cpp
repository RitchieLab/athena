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
#include "NormalizeContinuous.h"
#include <Descriptive.h>

#include <iostream>
using namespace std;

namespace data_manage
{

NormalizeContinuous::NormalizeContinuous()
{
}

NormalizeContinuous::~NormalizeContinuous()
{
}


///
/// Rescales the continuous variables by taking absolute difference
/// from the mean and dividing by the standard deviation
/// @param holder Dataholder with all dat
/// @param varIndex Variable to be scaled
///
void NormalizeContinuous::adjustContin(Dataholder* holder, unsigned int varIndex){
	stat::Descriptive calculator;

	vector<float> values(holder->numInds(), 0);

	unsigned int currInd;
	for(currInd = 0; currInd < holder->numInds(); currInd++){
		values[currInd] = holder->getInd(currInd)->getCovariate(varIndex);
	}

	float meanValue = calculator.mean(values);
	float standardDev = calculator.standard_dev(values);

	Individual* ind;
	// subtract each value from the mean and divide by the standard deviation
	for(currInd=0; currInd < holder->numInds(); currInd++){
		ind = holder->getInd(currInd);
		ind->setCovariate(varIndex, ((ind->getCovariate(varIndex) - meanValue) / standardDev));
	}

}



///
/// Rescales the continuous variables by taking absolute difference
/// from the mean and dividing by the standard deviation
/// @param holder Dataholder with all dat
///
void NormalizeContinuous::adjustStatus(Dataholder* holder){
	stat::Descriptive calculator;
	
	vector<float> values(holder->numInds(), 0);
	
	unsigned int currInd;
	for(currInd=0; currInd < holder->numInds(); currInd++){
		values[currInd] = holder->getInd(currInd)->getStatus();
	}
	meanVal = calculator.mean(values);
	stdDev = calculator.standard_dev(values);

	Individual* ind;
	// subtract each value from the mean and divide by the standard deviation
	for(currInd=0; currInd < holder->numInds(); currInd++){
		ind=holder->getInd(currInd);
		ind->setStatus((ind->getStatus()-meanVal)/stdDev);
	}
	
}

}
