#ifndef SCALEDATA_H_
#define SCALEDATA_H_

#include "Dataholder.h"
#include <string>

namespace data_manage
{

///
/// Abstract base class that performs scaling of data
///
class ScaleData
{
public:
  ScaleData();
  virtual ~ScaleData();

  virtual void adjust_contin(Dataholder* holder, unsigned int var_index){return;}
  
  virtual void adjust_status(Dataholder* holder){return;}
  
  virtual float get_original_status(float status){return status;}

  void sigmoid_status(Dataholder* holder);
  
  virtual std::string output_scale_info(){return "";}

};

}

#endif /*DATAADJUST_H_*/
