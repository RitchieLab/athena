#ifndef DESCRIPTIVE_H_
#define DESCRIPTIVE_H_

#include <cmath>
#include <vector>

namespace stat
{

class Descriptive
{
public:
  Descriptive();
  ~Descriptive();

  /// returns the mean for a vector of floats
  float mean(std::vector<float>& values);

  /// returns the standard deviation for a vector of floats
  float standard_dev(std::vector<float>& values);

};

}

#endif /*DESCRIPTIVE_H_*/
