#include "OttDummyConvert.h"

#include <iostream>
using namespace std;

namespace data_manage
{

OttDummyConvert::OttDummyConvert()
{
  convertor.assign(4, vector<char>(2,0));
  convertor[0][0] = -1;
  convertor[0][1] = -1;
  convertor[1][0] = 0;
  convertor[1][1] = 2;
  convertor[2][0] = 1;
  convertor[2][1] = -1;
  convertor[3][0] = 3;
  convertor[3][1] = 3;

}

OttDummyConvert::~OttDummyConvert()
{
}



///
/// Converts genotypes to ott dummy ones
///
void OttDummyConvert::convert_genotypes(Dataholder* holder){

  Individual* ind;
  vector<char> new_genos;
  unsigned int curr_loc;
  unsigned int num_loci = holder->num_genos();

  for(unsigned int curr_ind = 0; curr_ind < holder->num_inds(); curr_ind++){
    ind = holder->get_ind(curr_ind);

    new_genos.clear();
    for(curr_loc=0; curr_loc < num_loci; curr_loc++){
      // insert 2 loci that are recoded from original value
      new_genos.insert(new_genos.end(), convertor[ind->get_genotype(curr_loc)].begin(),
          convertor[ind->get_genotype(curr_loc)].end());
    }
    // replace the original genotypes with the new ones
    ind->set_all_genotypes(new_genos);
  }

  holder->ott_dummy_encoding(true);
}


}
