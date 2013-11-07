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
#include "Dataset.h"

using namespace std;

namespace data_manage
{

Dataset::Dataset()
{
	ssTotal = 0;
	binaryStatusOnly = false;
}

Dataset::~Dataset()
{
}


///
/// Appends pointer list to end of current list
/// @param newInds vector of Individual pointers
///
void Dataset::addInds(vector<Individual* >& newInds){
	inds.insert(inds.end(), newInds.begin(), newInds.end());
}


///
/// Outputs set displaying the status and genotypes for use in verifying
/// the sets
/// @param os ostream to write to
/// @param d Dataset containes current individuals
/// @return ostream
///
ostream& operator<<(ostream& os, Dataset& d){
	Individual * ind;
	for(unsigned int i=0; i<d.numInds(); i++){
		ind = d[i];
		os << ind->getID() << " ";
		os << ind->getStatus();
		for(unsigned int loc=0; loc < d.numGenos(); loc++){
			os << " " << ind->getGenotype(loc);
		}
		// output continuous after genotypes
		for(unsigned int cov=0; cov < d.numCovariates(); cov++){
			os << " " << ind->getCovariate(cov);
		}
		
		os << endl;
	}
	
	return os;
}


///
/// Calculates SSTotal on this set
///
void Dataset::calcSSTotal(){
	
	float diff=0.0;
	float meanVal;
	float statTotal = 0.0;
	
	for(unsigned int i=0; i < inds.size(); i++){
		statTotal += inds[i]->getStatus();
	}
	
	meanVal = statTotal/inds.size();
	
	for(unsigned int i=0; i<inds.size(); i++){
		diff = diff + (inds[i]->getStatus() - meanVal) * (inds[i]->getStatus() - meanVal);
	}
	
	ssTotal = diff;
}


///
/// Combines two datasets into one new one
///
Dataset Dataset::operator+(Dataset& d){
	Dataset newSet;
	newSet = *this;

	for(unsigned int i=0; i<d.numInds(); i++){
		newSet.addInd(d.getInd(i));
	}

	newSet.calcSSTotal();
	return newSet;
}


}
