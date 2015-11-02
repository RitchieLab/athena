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
#include "ScaleContinuous.h"

#include <sstream>
#include <iostream>
using namespace std;

namespace data_manage
{

ScaleContinuous::ScaleContinuous()
{
}

ScaleContinuous::~ScaleContinuous()
{
}


///
/// For all continuous variables, the values are all
/// scaled by dividing by the largest value in the variable
/// @param holder Dataholder with all dat
///
void ScaleContinuous::adjustContin(Dataholder* holder){
	for(unsigned int c=0; c < holder->numCovariates(); c++){
		adjustContin(holder, c);
	}
}

///
/// For the continuous variable in question, the values are all
/// scaled by dividing by the largest value in the variable
/// @param holder Dataholder with all dat
/// @param varIndex Variable to be scaled
///
void ScaleContinuous::adjustContin(Dataholder* holder, unsigned int varIndex){
	unsigned int currInd;
	Individual* ind;

	float maxValue = -1e30;
	float covarMin = 1e30;
	
	for(currInd = 0; currInd < holder->numInds(); currInd++){
		if(holder->getInd(currInd)->getCovariate(varIndex) == holder->getMissingCoValue()){
			continue;
		}
		if(holder->getInd(currInd)->getCovariate(varIndex) > maxValue)
			maxValue = holder->getInd(currInd)->getCovariate(varIndex);
		if(holder->getInd(currInd)->getCovariate(varIndex) < covarMin)
			covarMin = holder->getInd(currInd)->getCovariate(varIndex);
	}

	float covarDiff = maxValue-covarMin;

	// scale from 0 to 1
	for(currInd=0; currInd < holder->numInds(); currInd++){
		if(holder->getInd(currInd)->getCovariate(varIndex) == holder->getMissingCoValue())
			continue;
		ind = holder->getInd(currInd);
		// scaledValue = (rawValue - min) / (max - min);
		ind->setCovariate(varIndex, (ind->getCovariate(varIndex)-covarMin)/covarDiff);
	}

}


///
/// Divides all status values by the maximum value so they
/// will be scaled from zero to one.
/// @param holder Dataholder with all data
///
void ScaleContinuous::adjustStatus(Dataholder* holder){
	statMax = holder->getInd(0)->getStatus();
	statMin = holder->getInd(0)->getStatus();
	
	unsigned int currInd;
	Individual* ind;
	
	float status;
	
	for(currInd=0; currInd < holder->numInds(); currInd++){
		status = holder->getInd(currInd)->getStatus();
		if(status > statMax)
			statMax = status;
		if(status < statMin)
			statMin = status;
	}

	float statDiff = statMax-statMin;
	
	// divide all values by max and set status to that
	for(currInd=0; currInd < holder->numInds(); currInd++){
		ind=holder->getInd(currInd);
		ind->setStatus((ind->getStatus()-statMin)/statDiff);
	}  
}


///
/// Returns string that gives information on scaling performed
/// @return string
///
string ScaleContinuous::outputScaleInfo(){
	
	stringstream ss;
	ss << "ScaleMax=" << statMax << " StatusMin=" << statMin << std::endl;
	return ss.str();
}

}
