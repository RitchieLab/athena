#ifndef SAMPLENOREPLACEMENT_H_
#define SAMPLENOREPLACEMENT_H_

#include <cstdlib>
#include <vector>

namespace data_manage
{

///
/// Perform sampling without replacement
///

class RandomNoReplace
{
public:

  /// Sampled indexes are returned in the samples vector
  static void SampleWithoutReplacement(int popSize, int sampSize, std::vector<int>& samples);
  
private:
};

}

#endif /*SAMPLENOREPLACEMENT_H_*/
