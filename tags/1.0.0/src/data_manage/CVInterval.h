#ifndef CVINTERVAL_H_
#define CVINTERVAL_H_

#include "Dataset.h"

namespace data_manage
{

///
/// A crossvalidation interval with a dataset for a training, testing, and
/// an optional validation set.
///

class CVInterval
{
public:
  CVInterval();
  CVInterval(unsigned int num);
  ~CVInterval();

  Dataset& get_set(unsigned int index){return sets[index];}

  unsigned int num_sets(){return sets.size();}

  void num_sets(unsigned int num);

  void add_set(Dataset& set){sets.push_back(set);}

  Dataset& get_training(){return sets[0];}
  Dataset& get_testing(){return sets[1];}
  Dataset& get_validation(){return sets[2];}

private:

  std::vector<Dataset> sets;

};

}

#endif /*CVINTERVAL_H_*/
