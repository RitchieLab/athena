#ifndef SCALECONTINUOUS_H_
#define SCALECONTINUOUS_H_

#include "ScaleData.h"

namespace data_manage
{

class ScaleContinuous : public data_manage::ScaleData
{
public:
  ScaleContinuous();
  ~ScaleContinuous();

  void adjust_contin(Dataholder* holder, unsigned int var_index);
  
  /// Scales status by dividing all 
  void adjust_status(Dataholder* holder);

  /// Returns original value for a scaled status value
  float get_original_status(float status){return (status*statmax)-statadjust;}
  
  /// Writes scaling information to ostream
  string output_scale_info();
  
private:
  
  float statmin, statmax, statadjust;
  

};

}

#endif /*SCALECONTINUOUS_H_*/
