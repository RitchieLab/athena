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
#include "ScaledDataFactory.h"
#include "NormalizeContinuous.h"
#include "ScaleContinuous.h"

namespace data_manage{

map<string, ScaledDataFactory::ScaleType> ScaledDataFactory::scaleMap;


///
/// Function creates a data scaler based on the name passed
/// @param scaleType name of scaling type
/// @return pointer to new ScaleData
/// @throws DataExcept when no match
///
ScaleData* ScaledDataFactory::createScaler(string scaleType){
	
	if(scaleMap.empty()){
		setScaleMap();
	}
	
	ScaleData* newScaler;
	
	switch(scaleMap[scaleType]){
		case NoMatch:
			throw DataExcept("No scaler matching " + scaleType);
			break;
		case ScaleNorm:
			newScaler = new NormalizeContinuous;
			break;
		case ScaleContin:
			newScaler = new ScaleContinuous;
			break;
		case NoScale:
			newScaler = new ScaleData;
			break;
		default:
			throw DataExcept("No scaler matching " + scaleType);
	}
	
	return newScaler;
	
}


///
/// Establishes the map for use in creating solutions
/// @return 
///
void ScaledDataFactory::setScaleMap(){
	scaleMap["NORMSTDEV"]=ScaleNorm;
	scaleMap["NORMMAX"]=ScaleContin;
	scaleMap["NONE"]=NoScale;

}

}
