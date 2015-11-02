#include "Descriptive.h"

using namespace std;

namespace stat
{

Descriptive::Descriptive()
{
}

Descriptive::~Descriptive()
{
}

///
/// returns the mean for a vector of floats
/// @param values
/// @return mean
///
float Descriptive::mean(vector<float>& values){
  float mean_val=0.0;
  vector<float>::iterator iter;
  for(iter = values.begin(); iter != values.end(); ++iter)
    mean_val += *iter;

  return mean_val / values.size();

}


///
/// returns the standard deviation for a vector of floats
/// @param values
/// @return
///
float Descriptive::standard_dev(vector<float>& values){
  float mean_val = mean(values);

  vector<float>::iterator iter;
  float total=0.0;
  for(iter = values.begin(); iter != values.end(); ++iter){
    total += pow(fabs(*iter - mean_val), 2);
  }

  total = total / values.size();

  return sqrt(total);


}


}
