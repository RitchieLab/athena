#include "StephenDummyConvert.h"

#include <iostream>
using namespace std;

namespace data_manage
{

StephenDummyConvert::StephenDummyConvert()
{
}

StephenDummyConvert::~StephenDummyConvert()
{
}


///
/// Converts genotypes to ott dummy ones
///
void StephenDummyConvert::convert_genotypes(Dataholder* holder){

  Individual* ind;

  unsigned int curr_loc;
  unsigned int num_loci = holder->num_genos();

  for(unsigned int curr_ind = 0; curr_ind < holder->num_inds(); curr_ind++){

    ind = holder->get_ind(curr_ind);

    for(curr_loc=0; curr_loc < num_loci; curr_loc++){
      if(ind->get_genotype(curr_loc) < 3){
        ind->set_genotype(curr_loc, ind->get_genotype(curr_loc)-1);
      }

    }
  }

  // when it is ott dummy encoded there are 2 variables for each genotype
  // in this case we still only use one variable
  holder->ott_dummy_encoding(false);
}


}
