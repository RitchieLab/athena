#include "ContinFileReader.h"
#include "Stringmanip.h"
#include <fstream>
#include <sstream>

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
      
      while(!ss.eof()){
        ss >> value;
        holder->get_ind_by_id(ind_id)->add_covariate(value);
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
