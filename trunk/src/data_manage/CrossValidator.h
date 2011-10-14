#ifndef CROSSVALIDATOR_H_
#define CROSSVALIDATOR_H_

#include "CVSet.h"
#include "Dataholder.h"

namespace data_manage
{


class TRandom
{
  public:
  int operator()(int n)
  {
    return rand() % n;
  }

};

///
///  Performs crossvalidation on a dataholder to return datasets with the
///

class CrossValidator
{
public:
  CrossValidator();
  ~CrossValidator();

  /// splits data into testing and training sets for the indicated number of intervals
  CVSet split_data(unsigned int num_crossval, Dataholder* holder);

private:

  void distribute_inds(unsigned int num_splits, std::vector<Individual*>& inds,
      std::vector<std::vector<Individual*> >& splits);

  void shuffle_inds(std::vector<Individual*> & inds);

  void status_bin(Dataholder* holder, std::vector<Individual*>& affected,
      std::vector<Individual*>& unaffected);
      
  CVSet split_by_num(Dataholder* holder);    

};

}

#endif /*CROSSVALIDATOR_H_*/
