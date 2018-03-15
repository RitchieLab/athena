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
#include "StephenDummyConvert.h"

#include <iostream>
using namespace std;

namespace data_manage
{

StephenDummyConvert::StephenDummyConvert()
{
}

StephenDummyConvert::~StephenDummyConvert()
{
}


///
/// Converts genotypes to ott dummy ones
///
void StephenDummyConvert::convertGenotypes(Dataholder* holder){

	Individual* ind;

	unsigned int currLoc;
	unsigned int numLoci = holder->numGenos();

	for(unsigned int currInd = 0; currInd < holder->numInds(); currInd++){
		ind = holder->getInd(currInd);

		for(currLoc=0; currLoc < numLoci; currLoc++){
			if(ind->getGenotype(currLoc) < 3){
				ind->setGenotype(currLoc, ind->getGenotype(currLoc)-1);
			}

		}
	}

	// when it is ott dummy encoded there are 2 variables for each genotype
	// in this case we still only use one variable
	holder->ottDummyEncoding(false);
}

}
