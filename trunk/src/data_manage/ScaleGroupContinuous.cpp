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
#include "ScaleGroupContinuous.h"

#include <sstream>
#include <iostream>
using namespace std;

namespace data_manage
{

ScaleGroupContinuous::ScaleGroupContinuous()
{
}

ScaleGroupContinuous::~ScaleGroupContinuous()
{
}


///
/// For the continuous variable in question, the values are all
/// scaled by dividing by the largest value int the previously set group
/// @param holder Dataholder with all dat
/// @param varIndex Variable to be scaled
///
void ScaleGroupContinuous::adjustContin(Dataholder* holder, unsigned int varIndex){
	unsigned int currInd;
	Individual* ind;
	
	string groupName = holder->getCovarGroupName(varIndex);
	
	// divide all values by largest and set in dataholder
	for(currInd=0; currInd < holder->numInds(); currInd++){
		if(holder->getInd(currInd)->getCovariate(varIndex) == holder->getMissingCoValue())
			continue;
		ind = holder->getInd(currInd);
		ind->setCovariate(varIndex, (ind->getCovariate(varIndex)-minGroupValues[groupName] )/groupCovarDiff[groupName]);
	}

}


///
/// Determines maximum value for each continuous input group for use in scaling
/// inputs.
/// @param dataholder Dataholder
///
void ScaleGroupContinuous::adjustContin(Dataholder* holder){
	std::map<std::string, std::vector<int> > groupMap = holder->getContinGroups();
	std::map<std::string, std::vector<int> >::iterator groupIter;
	std::vector<int>::iterator continIter;
	unsigned int totalInds = holder->numInds();
	
	float missingCoVal = holder->getMissingCoValue();
	float continValue;
	
	// set values to absurdly small size
	maxGroupValues.clear();
	minGroupValues.clear();
	
	for(groupIter=groupMap.begin(); groupIter != groupMap.end(); ++groupIter){
		maxGroupValues[groupIter->first]=-1e30;
		minGroupValues[groupIter->first]=1e30;
		string groupName = groupIter->first;
		for(unsigned int ind=0; ind < totalInds; ind++){
			for(continIter=groupIter->second.begin(); continIter != groupIter->second.end(); continIter++){
				continValue = holder->getInd(ind)->getCovariate(*continIter);
				if(continValue != missingCoVal && continValue){
					if(continValue > maxGroupValues[groupName]){
						maxGroupValues[groupName] = continValue;
					}
					if(continValue < minGroupValues[groupName]){
						minGroupValues[groupName] = continValue;
					}
				}
			}
		}
	}
	
	groupCovarDiff.clear();
	
	// iterate through again and reset all continuous values
	for(groupIter=groupMap.begin(); groupIter != groupMap.end(); ++groupIter){
		string groupName = groupIter->first;
		groupCovarDiff[groupName] = maxGroupValues[groupName]-minGroupValues[groupName];
	
		for(unsigned int ind=0; ind < totalInds; ind++){
			for(continIter=groupIter->second.begin(); continIter != groupIter->second.end(); continIter++){
				continValue = holder->getInd(ind)->getCovariate(*continIter);
				if(continValue != missingCoVal){
					// scaledValue = (rawValue - min) / (max - min);
					holder->getInd(ind)->setCovariate(*continIter, 
						(continValue-minGroupValues[groupName]) / groupCovarDiff[groupName]);
				}
			}
		}
	}	
	
}


///
/// Divides all status values by the maximum value so they
/// will be scaled from zero to one.
/// @param holder Dataholder with all data
///
void ScaleGroupContinuous::adjustStatus(Dataholder* holder){
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
string ScaleGroupContinuous::outputScaleInfo(){
	
	stringstream ss;
	
	ss << "ScaleMax=" << statMax << " StatusMin==" << statMin << std::endl;
	return ss.str();
}

}
