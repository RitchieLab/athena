#ifndef OUTPUTSET_H_
#define OUTPUTSET_H_

#include "Dataset.h"
#include <string>

namespace data_manage
{

///
/// Outputs Dataset to file
///
class OutputSet
{

public:

  /// outputs designated set
  void outputSet(std::string name, Dataset& set);
  
  /// outputs set and adds cv number to name
  void outputCV(std::string name, Dataset& set, int cv);
  
private:

};

}

#endif /*DATASET_H_*/
