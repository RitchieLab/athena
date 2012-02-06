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
