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

#include "MDRFileHandler.h"
#include <fstream>
#include <sstream>
#include "Stringmanip.h"
#include <math.h>
#include <float.h>

#include <iostream>
#include <cstdlib>


namespace data_manage
{

MDRFileHandler::MDRFileHandler()
{
  dummy_id = 1;
}

MDRFileHandler::~MDRFileHandler()
{
}



///
///Parses the data file and fills this object
///@param datafile name of data file
///@param holder DataHolder to fill with
///data from the datafile
///@param missingValue value of missing data in current set
///@param statusMissingValue value that indicates missing status for an individual
///@param containds_id true when file has first column as ID value
///@return none
///@throws DataExcept on error
///
void MDRFileHandler::parse_file(string filename, Dataholder * holder, 
        int missingValue, float statusMissingValue, bool contains_id){
        
  int max_locus_value = 0;
  ifstream data_stream(filename.c_str(), ios::in);

  if(!data_stream.is_open()){
    throw DataExcept("Error:  Unable to open " + filename + "\n");
  }

  string line, out, ind_id;
  string::size_type last_number_pos;
  double ind_status;
  int geno;

  getline(data_stream, line);

  // format is staus followed by genotypes separated by white space
  bool any_missing = false;

  do{
    if(line.find_first_of("0123456789") == string::npos || line.find("#")==0){
      getline(data_stream, line);
      continue;
    }
    // remove windows carriage return
    last_number_pos = line.find_last_of("0123456789");
    line = line.substr(0,last_number_pos+1);

    stringstream ss(line);
    Individual ind;
    if(contains_id){
        ss >> ind_id;
    }
    else{
        ind_id = Stringmanip::itos(dummy_id++);
    }
    ind.set_id(ind_id);
    ss >> ind_status;
    
    // check to see if ind should be skipped because status is missing
    if(fabs(ind_status - statusMissingValue) > DBL_EPSILON){
    
      ind.set_status(ind_status);
      while(ss >> geno){
        if(geno == missingValue){
          ind.add_genotype(3);
          any_missing = true;
        }
        else if(geno > 2){
        throw DataExcept("All genotypes must be less than 3 or equal the missing value " +
            Stringmanip::itos(missingValue) + "\nIndividual " + Stringmanip::itos(dummy_id++) +
            " has the genotype " + Stringmanip::itos(geno));
        }
        else{
          if(geno > max_locus_value){
            max_locus_value = geno;
          }
          ind.add_genotype(geno);
        }
      }
      holder->add_ind(ind);
    }
    else{
	cout << "Skipping individual " << ind_id << " with status " << ind_status << " (STATUSMISSINGVALUE is set to " << statusMissingValue << ")" << endl;
    }

    getline(data_stream, line);
  }while(!data_stream.eof());

  data_stream.close();

  holder->set_max_locus_value(max_locus_value);
  holder->any_missing_genos(any_missing);
  holder->set_missing_genotype(3);

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
void MDRFileHandler::parse_file(std::string train_file, std::string test_file, Dataholder* holder,
	  int missingValue, float statusMissingValue, bool contains_id){
  
  parse_file(train_file, holder, missingValue, statusMissingValue, contains_id);
  // now the number of inds in the set will specify the split point for CV generation
  holder->set_test_split(holder->num_inds());
  
  // parse the testing file and add to dataholder
  parse_file(test_file, holder, missingValue, statusMissingValue, contains_id);
  
}



///
///Writes the data file to a text file using
///standard MDR format text file format.
///@param datafile name of data file
///@param dataholder DataHolder containing
///data to write to file
///@return none
///
void MDRFileHandler::write_file(string filename,
  Dataholder* dataholder){

  std::ofstream outStream(filename.c_str(), ios::out);
  if(!outStream.is_open()){
    throw DataExcept("ERROR: Unable to open " + filename + " for writing\n");
  }

  unsigned int numInds = dataholder->num_inds();
  unsigned int lastLocus = dataholder->num_genos();
  unsigned int currLoc;

  Individual* curr_ind;

  for(unsigned int currInd=0; currInd < numInds; currInd++){
    curr_ind = dataholder->get_ind(currInd);
    outStream << Stringmanip::itos(curr_ind->get_status())
      << " ";
    for(currLoc=0; currLoc < lastLocus; currLoc++){
      outStream << Stringmanip::itos(curr_ind->get_genotype(currLoc))
        << " ";
    }
    outStream << Stringmanip::itos(curr_ind->get_genotype(currLoc)) << endl;
  }

  outStream.close();
}





}
