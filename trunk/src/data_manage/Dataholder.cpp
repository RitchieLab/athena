#include "Dataholder.h"
#include "Stringmanip.h"

namespace data_manage
{


Dataholder::Dataholder()
{
  any_missing = false;
  ott_encoded = false;
  max_locus = 0;
  split_num = -1;
}


///
/// Destructor frees all memory
///
Dataholder::~Dataholder()
{
  vector<Individual*>::iterator iter;
  for(iter = inds.begin(); iter != inds.end(); iter++)
    delete *iter;
}



///
/// Adds individual to set
/// @param ind Individual to add
///
void Dataholder::add_ind(Individual& ind){
  Individual* new_ind = new Individual(ind);
  inds.push_back(new_ind);
  inds_map[new_ind->get_id()] = new_ind;
}


///
/// Adds individual to set
/// @param ind Individual to add
///
void Dataholder::add_ind(Individual* ind){
  inds.push_back(ind);
  inds_map[ind->get_id()] = ind;
}

///
/// Adds default snp names to set
///
void Dataholder::add_default_snps(){
   unsigned int total = num_genos();
   for(unsigned int i=1; i<=total; i++){
     add_geno_name(Stringmanip::itos(i));
   }
}

///
/// Adds default covariate names to holder
///
void Dataholder::add_default_covars(){
    unsigned int total = inds[0]->num_covariates();
    for(unsigned int i=1; i<=total; i++){
      add_covar_name(Stringmanip::itos(i));  
    }
}


}
