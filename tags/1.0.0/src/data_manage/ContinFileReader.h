#ifndef CONTINFILEREADER_H_
#define CONTINFILEREADER_H_

#include "Dataholder.h"
#include "DataExcept.h"

namespace data_manage
{

///
/// Reads continuous variables from file and stores in Dataholder.  The
/// order of individuals must match that of the original data file.
///
class ContinFileReader
{
public:
  ContinFileReader();
  ~ContinFileReader();

  void read_contin_file(std::string filename, Dataholder* holder, float missingValue,
          bool contains_id=false);
          
  void read_contin_file(std::string train_file, std::string test_file, Dataholder* holder, 
    float missingValue, bool contains_id=false);  

private:
  int dummy_id;
};

}

#endif /*CONTINFILEREADER_H_*/
