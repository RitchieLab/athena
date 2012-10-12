/*
Copyright Marylyn Ritchie 2011

This file is part of ATHENA.

ATHENA is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

ATHENA is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with ATHENA.  If not, see <http://www.gnu.org/licenses/>.
*/
#include "TransferData.h"

#ifdef PARALLEL

using namespace data_manage;

///
/// sends all data to slaves 
/// @param holder Dataholder that contains data to send
///
void TransferData::sendData(data_manage::Dataholder& holder){

  // send each individual in dataset in original order
  unsigned int num_inds = holder.num_inds();
  
  Individual* ind;
  
  // broadcast number of inds, number of genotypes, and number of covariates
  // to all slaves
  unsigned int num_genotypes = holder.num_genos();
  unsigned int num_covariates = holder.num_covariates();
  
  int * set_info = new int[8];
  set_info[0] = num_inds;
  set_info[1] = num_genotypes;
  set_info[2] = num_covariates;
  set_info[3] = holder.get_missing_genotype();
  set_info[4] = holder.ott_dummy_encoding();
  set_info[5] = holder.any_missing_genos();
  set_info[6] = holder.get_max_locus_value();
  set_info[7] = holder.get_test_split();
  
  MPI_Bcast(set_info, 7, MPI_INT, 0, MPI_COMM_WORLD);
  
  delete [] set_info;
  
  float missing_covalue;
  missing_covalue = holder.get_missing_covalue();

  MPI_Bcast(&missing_covalue, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
  holder.set_missing_covalue(missing_covalue);
  
  int id_string_size = 40;
  
  // establish array to hold covariates/status
  float * costatus = new float[num_covariates+1];
  // establish array to hold genotypes for sending
  char * genotypes = new char[num_genotypes+id_string_size];
  
  for(unsigned int i=0; i<num_inds; i++){
    ind = holder[i];
    // send status and all covariates 
    costatus[0] = ind->get_status();
    for(unsigned int covar=0; covar < num_covariates; covar++){
      costatus[covar+1] = ind->get_covariate(covar);
    }
    
    // send float values for this individual
    MPI_Bcast(costatus, num_covariates+1, MPI_FLOAT, 0, MPI_COMM_WORLD);
    // fill genotypes array and send along with the id string at end
    unsigned int geno;
    for(geno=0; geno < num_genotypes; geno++){
      genotypes[geno] = ind->get_genotype(geno);
    }
    
    string idstring = ind->get_id();
    // add string ID to genotypes list
    for(unsigned int id=0; id < idstring.size(); id++){
      genotypes[geno++] = idstring[id];
    }
    genotypes[geno] = '\0';
    
    MPI_Bcast(genotypes, num_genotypes, MPI_CHAR, 0, MPI_COMM_WORLD);
  }
  
  delete [] costatus;
  delete [] genotypes;
}
  

///
/// receives all data from slaves
/// @param holder Dataholder that will contain data
///
void TransferData::receiveData(data_manage::Dataholder& holder){

// cout << "receiveData" << endl;
  
  int * set_info = new int[8];
  
  MPI_Bcast(set_info, 8, MPI_INT, 0, MPI_COMM_WORLD);
  
  int num_inds = set_info[0];
  int num_genotypes = set_info[1];
  int num_covariates = set_info[2];
  holder.set_missing_genotype(set_info[3]);
  holder.ott_dummy_encoding(set_info[4]);
  holder.any_missing_genos(set_info[5]);
  holder.set_max_locus_value(set_info[6]);
  holder.set_test_split(set_info[7]);

// cout << "Slave set_info: " << num_inds << " " << num_genotypes << " " << num_covariates <<
//   " " << holder.get_missing_genotype() << " " << holder.ott_dummy_encoding() << " "
//   << holder.any_missing_genos() <<  " " << holder.get_max_locus_value() << endl;

  delete [] set_info;
  
  float missing_covalue;
  MPI_Bcast(&missing_covalue, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
  holder.set_missing_covalue(missing_covalue);

// cout << "Slave missing covalue=" << holder.get_missing_covalue() << endl;
  
  int id_string_size = 40;
  
  // establish array to hold covariates/status
  float * costatus = new float[num_covariates+1];
  // establish array to hold genotypes for sending
  char * genotypes = new char[num_genotypes+id_string_size];  
  
  Individual* ind;
  int geno;
  
  for(int i=0; i<num_inds; i++){
    ind = new Individual;
    // receive float values for this individual
    MPI_Bcast(costatus, num_covariates+1, MPI_FLOAT, 0, MPI_COMM_WORLD);

    ind->set_status(costatus[0]);
    for(int covar=1; covar <= num_covariates; covar++){
      ind->add_covariate(costatus[covar]);
    }

    // receive genotypes for this individual
    MPI_Bcast(genotypes, num_genotypes, MPI_CHAR, 0, MPI_COMM_WORLD);
    // add genotypes
    for(geno=0; geno < num_genotypes; geno++){
      ind->add_genotype(genotypes[geno]);
    }
 
    // add the id string 
    string idstring;
    for(unsigned int id=0; id < idstring.size(); id++){
      if(genotypes[geno] == '\0')
        break;
      idstring += genotypes[geno++];
    }
    ind->set_id(idstring);
    holder.add_ind(ind);
  }
  
  delete [] costatus;
  delete [] genotypes;

}

#endif /* end if PARALLEL defined */

