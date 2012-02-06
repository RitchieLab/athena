#include "CVInterval.h"

namespace data_manage
{

///
/// Constructor
///
CVInterval::CVInterval()
{
}


///
/// Alternative constructor
/// @param num Number of sets that the individuals will be split into
///
CVInterval::CVInterval(unsigned int num){
  num_sets(num);
}

///
/// Destructor
///
CVInterval::~CVInterval()
{
}


///
/// Sets number of sets
/// @param num Number of sets that the individuals will be split into
///
void CVInterval::num_sets(unsigned int num){

  Dataset blank;
  /// fills with blank sets
  sets.assign(num, blank);
}




}
