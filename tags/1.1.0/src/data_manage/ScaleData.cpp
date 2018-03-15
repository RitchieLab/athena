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
#include "ScaleData.h"
#include <cmath>

namespace data_manage
{

ScaleData::ScaleData()
{
}

ScaleData::~ScaleData()
{
}

///
/// Alters status by putting status through a sigmoid function
/// @param holder Dataholder
///
void ScaleData::sigmoidStatus(Dataholder* holder){
	
	Individual * ind; 
	float status;
	
	for(unsigned int currInd = 0; currInd < holder->numInds(); currInd++){
		ind = (*holder)[currInd];
		status = ind->getStatus();
		if(status < -709) ind->setStatus(-1.0);
		else if(status > 709)ind->setStatus(1.0);
		else ind->setStatus(1.0/(1.0+exp(-status)));
	}
				
}


}
