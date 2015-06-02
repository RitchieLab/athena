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
#include "ScaleCategorical.h"

#include <sstream>
#include <iostream>
#include <set>
using namespace std;

namespace data_manage
{

ScaleCategorical::ScaleCategorical()
{
}

ScaleCategorical::~ScaleCategorical()
{
}


///
/// Values are converted into categorical values starting with 0.
/// @param holder Dataholder with all dat
///
void ScaleCategorical::adjustContin(Dataholder* holder){
	for(unsigned int c=0; c < holder->numCovariates(); c++){
		adjustContin(holder, c);
	}
}

///
/// Values are converted into categorical values starting with 0.
/// @param holder Dataholder with all dat
///
void ScaleCategorical::adjustContin(Dataset* dataset){
	for(unsigned int c=0; c < dataset->numCovariates(); c++){
		adjustContin(dataset, c);
	}
}


///
/// For the variable specified the number of different values are checked and
/// then the values are re-done starting with 0
/// @param holder Dataholder with all dat
/// @param varIndex Variable to be scaled
///
void ScaleCategorical::adjustContin(Dataholder* holder, unsigned int varIndex){
	unsigned int currInd, numInds = holder->numInds();
	Individual* ind;

	set<int> uniqueValues;

	for(currInd = 0; currInd < numInds; currInd++){
		ind = holder->getInd(currInd);
		if(ind->getCovariate(varIndex) == holder->getMissingCoValue()){
			continue;
		}
		int value = int(ind->getCovariate(varIndex));
		if(uniqueValues.find(value) == uniqueValues.end()){
			uniqueValues.insert(value);
		}
	}
	
	map<int, int> conversion;
	int conv=0;
	for(set<int>::iterator iter=uniqueValues.begin(); iter != uniqueValues.end();
		++iter){
		conversion[*iter]=conv;
		conv++;
	}

	for(currInd = 0; currInd < numInds; currInd++){
		ind = holder->getInd(currInd);
		if(ind->getCovariate(varIndex) == holder->getMissingCoValue()){
			continue;
		}
		int oldvalue = ind->getCovariate(varIndex);
		ind->setCovariate(varIndex, conversion[oldvalue]);		
	}
	
// 	holder->setNumLevels(varIndex, conversion.size());
	
}


///
/// For the variable specified the number of different values are checked and
/// then the values are re-done starting with 0
/// @param holder Dataset with all dat
/// @param varIndex Variable to be scaled
///
void ScaleCategorical::adjustContin(Dataset* dataset, unsigned int varIndex){
	unsigned int currInd, numInds = dataset->numInds();
	Individual* ind;
	set<int> uniqueValues;

	for(currInd = 0; currInd < numInds; currInd++){
		ind = dataset->getInd(currInd);
		if(ind->getCovariate(varIndex) == dataset->getMissingCoValue()){
			continue;
		}
		int value = int(ind->getCovariate(varIndex));
		if(uniqueValues.find(value) == uniqueValues.end()){
			uniqueValues.insert(value);
		}
	}
	
	map<int, int> conversion;
	int conv=0;
	for(set<int>::iterator iter=uniqueValues.begin(); iter != uniqueValues.end();
		++iter){
		conversion[*iter]=conv;
		conv++;
	}
	for(currInd = 0; currInd < numInds; currInd++){
		ind = dataset->getInd(currInd);
		if(ind->getCovariate(varIndex) == dataset->getMissingCoValue()){
			continue;
		}
		int oldvalue = ind->getCovariate(varIndex);
		ind->setCovariate(varIndex, conversion[oldvalue]);
	}
	dataset->setNumLevels(varIndex, conversion.size());
}


///
/// Divides all status values by the maximum value so they
/// will be scaled from zero to one.
/// @param holder Dataholder with all data
///
void ScaleCategorical::adjustStatus(Dataholder* holder){
	unsigned int currInd, numInds = holder->numInds();
	Individual* ind;
	set<int> uniqueValues;
	originalVals.clear();
	
	for(currInd = 0; currInd < numInds; currInd++){
		ind = holder->getInd(currInd);
		int value = int(ind->getStatus());
		if(uniqueValues.find(value) == uniqueValues.end()){
			uniqueValues.insert(value);
			originalVals[value]=ind->getStatus();
		}
	}
	
	
	if(uniqueValues.size() ==2 && uniqueValues.find(1) != uniqueValues.end() && uniqueValues.find(0) != 
		uniqueValues.end()){
		return;
	}
	else if(uniqueValues.size()==1 && uniqueValues.find(1) != uniqueValues.end()){
		return;
	}
	
	map<int, int> conversion;
	int conv=0;
	for(set<int>::iterator iter=uniqueValues.begin(); iter != uniqueValues.end();
		++iter){
		conversion[*iter]=conv;
		conv++;
	}

	for(currInd = 0; currInd < numInds; currInd++){
		ind = holder->getInd(currInd);
		int oldStatus = ind->getStatus();
		ind->setStatus(conversion[oldStatus]);		
	}
	
// 	holder->setNumStatusLevels(conversion.size());
	
}


///
/// Scales data to be categorical
/// @param holder Dataset 
///
void ScaleCategorical::adjustStatus(Dataset* dataset){
	unsigned int currInd, numInds = dataset->numInds();
	Individual* ind;
	set<int> uniqueValues;
	originalVals.clear();
	for(currInd = 0; currInd < numInds; currInd++){
		ind = dataset->getInd(currInd);
		int value = int(ind->getStatus());
		if(uniqueValues.find(value) == uniqueValues.end()){
			uniqueValues.insert(value);
			originalVals[value]=ind->getStatus();
		}
	}
	
	if(uniqueValues.size() ==2 && uniqueValues.find(1) != uniqueValues.end() && uniqueValues.find(0) != 
		uniqueValues.end()){
		return;
	}
	else if(uniqueValues.size()==1 && uniqueValues.find(1) != uniqueValues.end()){
		return;
	}	
	
	map<int, int> conversion;
	int conv=0;
	for(set<int>::iterator iter=uniqueValues.begin(); iter != uniqueValues.end();
		++iter){
		conversion[*iter]=conv;
		conv++;
	}

	for(currInd = 0; currInd < numInds; currInd++){
		ind = dataset->getInd(currInd);
		int oldStatus = ind->getStatus();
		ind->setStatus(conversion[oldStatus]);		
	}
	
	dataset->setNumStatusLevels(conversion.size());
	
}


///
/// Returns string that gives information on scaling performed
/// @return string
///
string ScaleCategorical::outputScaleInfo(){
	
	stringstream ss;
	ss << "Converted to categorical" << endl;
	return ss.str();
}

}
