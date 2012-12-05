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
#include "ContinFileReader.h"
#include "Stringmanip.h"
#include <fstream>
#include <sstream>
#include <iostream>

namespace data_manage
{

ContinFileReader::ContinFileReader()
{
  dummy_id = 1;
}

ContinFileReader::~ContinFileReader()
{
}

///
/// Reads continuous variables from file and stores in Dataholder
/// @param filename string
/// @param holder Dataholder
/// @param missingValue value that identifies absence of data
/// @param contains_id true when first column is an ID that should match the
/// ID in the corresponding genotype data file
///
void ContinFileReader::read_contin_file(string filename, Dataholder* holder,
        float missingValue, bool contains_id){

    holder->set_missing_covalue(missingValue);

    std::ifstream c_stream(filename.c_str(), ios::in);

    if(!c_stream.is_open()){
      throw DataExcept("ERROR: Unable to open " + filename + "\n");
    }

    string line, ind_id;
    float value;
    unsigned int curr_ind=0;
    unsigned int num_covars=0;

    while(!c_stream.eof()){
      getline(c_stream, line);

      if(line.find_first_of("0123456789") == string::npos){
        continue;
      }

      stringstream ss(line);
      
      if(contains_id){
          ss >> ind_id;
      }
      else{
          ind_id = Stringmanip::itos(dummy_id++);
      }
     Individual * ind;
     try{
       ind = holder->get_ind_by_id(ind_id);
     }
     catch(DataExcept& de){
       cout << "Skipping individual " << ind_id << " while reading " << filename << endl;  
       continue;
     }
      
      while(ss >> value){
//        holder->get_ind_by_id(ind_id)->add_covariate(value);
        ind->add_covariate(value);
      }
      if(num_covars > 0 && num_covars != ind->num_covariates()){
        throw DataExcept("ERROR: in file " +  filename +" individual " + ind_id +" has " +Stringmanip::itos(ind->num_covariates()) + " values but previous line has " + Stringmanip::itos(num_covars));
      }
      else{
//         num_covars = holder->get_ind_by_id(ind_id)->num_covariates();
         num_covars = ind->num_covariates();
      }
      curr_ind++;
    }

    c_stream.close();
    
}


///
///Parses the data file and fills this object
///@param datafile name of data file
///@param holder DataHolder to fill with
///data from the datafile
///@param missingValue value of missing data in current set
///@param containds_id true when file has first column as ID value
///@return none
///@throws DataExcept on error
///
void ContinFileReader::read_contin_file(std::string train_file, std::string test_file, Dataholder* holder,
	  float missingValue, bool contains_id){
	  
  // parse training file and add to dataholder  
  read_contin_file(train_file, holder, missingValue, contains_id);
  
  // parse the testing file and add to dataholder
  read_contin_file(test_file, holder, missingValue, contains_id);
  
}

}
