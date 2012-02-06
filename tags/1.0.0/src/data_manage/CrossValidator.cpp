#include "CrossValidator.h"

#include <algorithm>
#include <stdlib.h>

namespace data_manage
{

CrossValidator::CrossValidator()
{
}

CrossValidator::~CrossValidator()
{
}


///
///  Splits data into testing and training sets for the number of
///  CV intervals selected
///  @param num_crossval Number of crossvalidation intervals
///  @param holder Dataholder containing data for this analysis
///	 @return CVSet
///
CVSet CrossValidator::split_data(unsigned int num_crossval, Dataholder* holder){

  if(holder->get_test_split() > 0)
    return split_by_num(holder);

  CVSet set;
  
  vector<Individual*> shuffled, affected, unaffected;
  status_bin(holder, affected, unaffected);
  shuffle_inds(affected);
  shuffle_inds(unaffected);
    
  vector<Individual*> temp;
  vector<vector<Individual*> > splits(num_crossval, temp);
  distribute_inds(num_crossval, affected, splits);
  distribute_inds(num_crossval, unaffected, splits);

  if(num_crossval > 1){
    unsigned int group;
    // using the splits stored in vector construct the CV Intervals and fill the set
    for(unsigned int curr_cv=0; curr_cv < num_crossval; curr_cv++){
      Dataset training(holder->get_missing_covalue(), holder->get_missing_genotype()), 
        testing(holder->get_missing_covalue(), holder->get_missing_genotype());
      for(group=0; group < num_crossval; group++){
        if(group != curr_cv)
          training.add_inds(splits[group]);
        else
          testing.add_inds(splits[group]);
      }
      CVInterval interval;
      training.calc_sstotal();
      testing.calc_sstotal();
      interval.add_set(training);
      interval.add_set(testing);
      set.add_interval(interval);
    }
  }
  else{ // only one interval so don't split data
    CVInterval interval;
    Dataset training(holder->get_missing_covalue(), holder->get_missing_genotype());
    training.add_inds(splits[0]);
    training.calc_sstotal();
    interval.add_set(training);
    set.add_interval(interval);
  }

  return set;
}


///
/// Splits data according to index in the holder class.
/// For use with user-specified splits.
/// @param holder Dataholder
///
CVSet CrossValidator::split_by_num(Dataholder* holder){
  int split_index = holder->get_test_split();
  CVSet set;

  Dataset training(holder->get_missing_covalue(), holder->get_missing_genotype()), 
    testing(holder->get_missing_covalue(), holder->get_missing_genotype());
  
  vector<Individual*> test_set, train_set;
  
  int i=0;
  for(i=0; i<split_index; i++){
    train_set.push_back(holder->get_ind(i));
  }  
  
  int nInds = int(holder->num_inds());
  for(; i<nInds; i++){
    test_set.push_back(holder->get_ind(i));
  }

  training.add_inds(train_set);
  testing.add_inds(test_set);
  
  CVInterval interval;
  training.calc_sstotal();
  testing.calc_sstotal();
  interval.add_set(training);
  interval.add_set(testing);
  set.add_interval(interval);
  
  return set;
}

///
/// Splits vector into specified set
/// @param num_splits
/// @param inds vector of Individual pointers
/// @param splits two-dimensional vector of Individual pointers that will contain
/// individual splits
///
void CrossValidator::distribute_inds(unsigned int num_splits, vector<Individual*>& inds,
    vector<vector<Individual*> >& splits){

  for(unsigned int i=0; i<inds.size(); i++){
    splits[i%num_splits].push_back(inds[i]);
  }

}


///
/// Distributes individuals between affected and unaffected arrays
/// @param holder Dataholder
/// @param affected Vector will contain pointers to affected inds
/// @param unaffected Vector will contain pointers to unaffected inds
///
void CrossValidator::status_bin(Dataholder* holder, vector<Individual*>& affected,
    vector<Individual*>& unaffected){

  affected.clear();
  unaffected.clear();

  unsigned int n_inds = holder->num_inds();
  Individual* ind;
  for(unsigned int i=0; i < n_inds; i++){
    ind = holder->get_ind(i);
    if(ind->get_status())
      affected.push_back(ind);
    else
      unaffected.push_back(ind);
  }

}



///
/// Shuffles individuals in the dataset
/// so that the splitting of the data will be independent of the
/// original order in the data file
/// @param indexes vector will contain pointers to now shuffled individuals
///
void CrossValidator::shuffle_inds(vector<Individual*> & inds){

  TRandom myrand;
  random_shuffle(inds.begin(), inds.end(), myrand);
}



}
