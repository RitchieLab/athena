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

map<string, ScaledDataFactory::ScaleType> ScaledDataFactory::ScaleMap;


///
/// Function creates a data scaler based on the name passed
/// @param scale_type name of scaling type
/// @return pointer to new ScaleData
/// @throws DataExcept when no match
///
ScaleData* ScaledDataFactory::create_scaler(string scale_type){
  
  if(ScaleMap.empty()){
    setScaleMap();
  }
  
  ScaleData* new_scaler;
  
  switch(ScaleMap[scale_type]){
    case NoMatch:
      throw DataExcept("No scaler matching " + scale_type);
      break;
    case ScaleNorm:
      new_scaler = new NormalizeContinuous;
      break;
    case ScaleContin:
      new_scaler = new ScaleContinuous;
      break;
    case NoScale:
      new_scaler = new ScaleData;
      break;
    default:
      throw DataExcept("No scaler matching " + scale_type);
  }
  
  return new_scaler;
  
}

///
/// Establishes the map for use in creating solutions
/// @return 
///
void ScaledDataFactory::setScaleMap(){
  ScaleMap["NORMSTDEV"]=ScaleNorm;
  ScaleMap["NORMMAX"]=ScaleContin;
  ScaleMap["NONE"]=NoScale;

}

}
