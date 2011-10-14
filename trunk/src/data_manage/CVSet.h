#ifndef CVSET_H_
#define CVSET_H_

#include "CVInterval.h"

namespace data_manage
{

class CVSet
{
public:
  CVSet();
  ~CVSet();

  /// Get number of intervals in set
  unsigned int num_intervals(){return intervals.size();}

  /// Get specific interval from set
  CVInterval&  get_interval(unsigned int index){return intervals[index];}

  /// Add interval to set
  void add_interval(CVInterval& interval){intervals.push_back(interval);}


private:
  std::vector<CVInterval> intervals;

};

}

#endif /*CVSET_H_*/
