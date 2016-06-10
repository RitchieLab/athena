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
#include "Dataholder.h"
#include "Stringmanip.h"
#include <Descriptive.h>
#include <iostream>

namespace data_manage
{


Dataholder::Dataholder()
{
	anyMissing = false;
	ottEncoded = false;
	binaryStatusOnly = false;
	maxLocus = 0;
	splitNum = -1;
	for(int i=0; i<10000; i++){
		statusConvert[i]=i;
	}
}


///
/// Destructor frees all memory
///
Dataholder::~Dataholder()
{
	vector<Individual*>::iterator iter;
	for(iter = inds.begin(); iter != inds.end(); iter++)
		delete *iter;
}



///
/// Adds individual to set
/// @param ind Individual to add
///
void Dataholder::addInd(Individual& ind){
	Individual* newInd = new Individual(ind);
	inds.push_back(newInd);
	indsMap[newInd->getID()] = newInd;
}


///
/// Adds individual to set
/// @param ind Individual to add
///
void Dataholder::addInd(Individual* ind){
	inds.push_back(ind);
	indsMap[ind->getID()] = ind;
}

///
/// Adds default snp names to set
///
void Dataholder::addDefaultSnps(){
	 unsigned int total = numGenos();
	 for(unsigned int i=1; i<=total; i++){
// 		 addGenoName(Stringmanip::itos(i));
		 addGenoName(Stringmanip::numberToString(i));
	 }
}

///
/// Adds default covariate names to holder
///
void Dataholder::addDefaultCovars(){
		unsigned int total = inds[0]->numCovariates();
		for(unsigned int i=1; i<=total; i++){
// 			addCovarName(Stringmanip::itos(i));
			addCovarName(Stringmanip::numberToString(i));
		}
}


///
/// Check for variance across all data splits
/// @param cvSet
///
void Dataholder::checkVariance(CVSet& cvSet){

	unsigned int nIntervals = cvSet.numIntervals();

	if(nIntervals == 1){
		// check entire set at one time
		checkVariance();
	}
	else{
		excludedGenos.clear();
		excludedContin.clear();
		for(unsigned int i=0; i<nIntervals; i++){
			CVInterval cvInt = cvSet.getInterval(i);
			Dataset trainSet = cvInt.getTraining();
			Dataset testSet = cvInt.getTesting();
			checkVariance(trainSet);
			checkVariance(testSet);
		}
	}

}


///
/// Check variance for the dataSet passed
///
void Dataholder::checkVariance(Dataset& dataSet){
	unsigned int nInds = dataSet.numInds();
	unsigned int nGenos = dataSet.numGenos();
	unsigned int nContin = dataSet.numCovariates();
	stat::Descriptive devCalculator;
	vector<float> vals(nInds, 0.0);

	// check each geno
	for(unsigned int i=0; i<nGenos; i++){
		for(unsigned int j=0; j<nInds; j++){
			vals[j] = dataSet[j]->getGenotype(i);
		}
		float stdDev = devCalculator.standard_dev(vals);
		if(stdDev == 0){
			excludedGenos.push_back(i);
		}
	}

	// check each continuous variable
	for(unsigned int i=0; i<nContin; i++){
		for(unsigned int j=0; j<nInds; j++){
			vals[j] = dataSet[j]->getCovariate(i);
		}
		float stdDev = devCalculator.standard_dev(vals);
		if(stdDev == 0){
			excludedContin.push_back(i);
		}
	}
}



///
/// Check for variance in variables and create lists of ones that have zero variance
///
void Dataholder::checkVariance(){
	excludedGenos.clear();
	excludedContin.clear();

	unsigned int nInds = numInds();
	unsigned int nGenos = numGenos();
	unsigned int nContin = numCovariates();
	stat::Descriptive devCalculator;
	vector<float> vals(nInds, 0.0);

	// check each geno
	for(unsigned int i=0; i<nGenos; i++){
		for(unsigned int j=0; j<nInds; j++){
			vals[j] = inds[j]->getGenotype(i);
		}
		float stdDev = devCalculator.standard_dev(vals);
		if(stdDev == 0){
			excludedGenos.push_back(i);
		}
	}

	// check each continuous variable
	for(unsigned int i=0; i<nContin; i++){
		for(unsigned int j=0; j<nInds; j++){
			vals[j] = inds[j]->getCovariate(i);
		}
		float stdDev = devCalculator.standard_dev(vals);
		if(stdDev == 0){
			excludedContin.push_back(i);
		}
	}
}

///
/// orders status as 0,1,2
/// saves map for converting back
///
void Dataholder::reNumberStatus(){
	map<int, int> statusMap;
	statusConvert.clear();
	int newValue=0;
	for(size_t i=0; i<inds.size(); i++){
		if(statusMap.find(int(inds[i]->getStatus())) == statusMap.end()){
			statusMap[int(inds[i]->getStatus())]=1;
		}
	}

	for(map<int, int>::iterator iter=statusMap.begin(); iter != statusMap.end();
		++iter){
		statusConvert[newValue]=iter->first;
		iter->second = newValue++;
	}

	for(size_t i=0; i<inds.size(); i++){
		inds[i]->setStatus(statusMap[int(inds[i]->getStatus())]);
	}
}

int Dataholder::getOriginalStatus(int status){
	return statusConvert[status];
}

///
/// Returns highest value in a continuous variable
/// @param varIndex index of variable to get value
///
float Dataholder::getHighestVal(int varIndex){
	return highContin[varIndex];
}


}
