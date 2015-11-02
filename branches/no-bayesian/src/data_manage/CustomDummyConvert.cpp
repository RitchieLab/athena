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
#include "CustomDummyConvert.h"
#include <algorithm>  

using namespace std;

namespace data_manage
{

CustomDummyConvert::CustomDummyConvert(vector<int>& tokens)
{
	for(vector<int>::iterator t=tokens.begin(); t!=tokens.end(); t++){
		convertor.push_back(*t);
	}
}

CustomDummyConvert::CustomDummyConvert()
{
}

CustomDummyConvert::~CustomDummyConvert()
{
}


///
/// Converts genotypes to ott dummy ones
///
void CustomDummyConvert::convertGenotypes(Dataholder* holder){

	Individual* ind;

	unsigned int currLoc;
	unsigned int numLoci = holder->numGenos();
	int missingGeno = holder->getMissingGenotype();
	int newMissing = missingGeno;

	// if the current missing genotype is part of conversion 
	// need to select a value that isn't part of it
	if(find(convertor.begin(), convertor.end(), missingGeno) != convertor.end()){
		newMissing = -255;
		while(find(convertor.begin(), convertor.end(), newMissing)!= convertor.end()){
			newMissing++;
		}
	}

	for(unsigned int currInd = 0; currInd < holder->numInds(); currInd++){
		ind = holder->getInd(currInd);

		for(currLoc=0; currLoc < numLoci; currLoc++){
			if(ind->getGenotype(currLoc) != missingGeno){
				ind->setGenotype(currLoc, convertor[ind->getGenotype(currLoc)]);
			}
			else{
				ind->setGenotype(currLoc, newMissing);
			}
		}
	}

	// in this case we still only use one variable
	holder->ottDummyEncoding(false);
	holder->setMissingGenotype(newMissing);
	
}


}
