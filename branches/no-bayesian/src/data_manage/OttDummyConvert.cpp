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
#include "OttDummyConvert.h"

#include <iostream>
using namespace std;

namespace data_manage
{

OttDummyConvert::OttDummyConvert()
{
	convertor.assign(4, vector<char>(2,0));
	convertor[0][0] = -1;
	convertor[0][1] = -1;
	convertor[1][0] = 0;
	convertor[1][1] = 2;
	convertor[2][0] = 1;
	convertor[2][1] = -1;
	convertor[3][0] = 3;
	convertor[3][1] = 3;

}

OttDummyConvert::~OttDummyConvert()
{
}



///
/// Converts genotypes to ott dummy ones
///
void OttDummyConvert::convertGenotypes(Dataholder* holder){
	Individual* ind;
	vector<char> newGenos;
	unsigned int currLoc;
	unsigned int numLoci = holder->numGenos();

	for(unsigned int currInd = 0; currInd < holder->numInds(); currInd++){
		ind = holder->getInd(currInd);

		newGenos.clear();
		for(currLoc=0; currLoc < numLoci; currLoc++){
			// insert 2 loci that are recoded from original value
			newGenos.insert(newGenos.end(), convertor[ind->getGenotype(currLoc)].begin(),
					convertor[ind->getGenotype(currLoc)].end());
		}
		// replace the original genotypes with the new ones
		ind->setAllGenotypes(newGenos);
	}

	setMultipleVars(true);
	holder->ottDummyEncoding(true);
}


}
