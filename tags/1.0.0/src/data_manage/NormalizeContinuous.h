#ifndef NORMALIZECONTINUOUS_H_
#define NORMALIZECONTINUOUS_H_

#include "ScaleData.h"


namespace data_manage
{

class NormalizeContinuous : public data_manage::ScaleData
{
public:
  NormalizeContinuous();
  ~NormalizeContinuous();

  void adjust_contin(Dataholder* holder, unsigned int var_index);
  
  /// Adjust continuous status values
  void adjust_status(Dataholder* holder);

  /// Returns original value for a scaled status value
  float get_original_status(float status){return (status*stddev)+meanval;}
  
private:
  
  float stddev, meanval;
  
};

}

#endif /*NORMALIZECONTINUOUS_H_*/
