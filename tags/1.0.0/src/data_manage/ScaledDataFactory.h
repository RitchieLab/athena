#ifndef SCALEDATAFACTORY_H_
#define SCALEDATAFACTORY_H_

#include "ScaleData.h"
#include "DataExcept.h"

namespace data_manage{

class ScaledDataFactory{

  public:

    static ScaleData* create_scaler(string scale_type);
    
  private:
    
    enum ScaleType{
      NoMatch,
      ScaleNorm,
      ScaleContin,
      NoScale
    };

    static map<string, ScaleType> ScaleMap;
    
    /// Sets the map
    static void setScaleMap();
    
};


}

#endif /*SCALEDATAFACTORY_H_*/
