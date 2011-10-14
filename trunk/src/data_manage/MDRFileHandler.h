#ifndef MDRFILEHANDLER_H_
#define MDRFILEHANDLER_H_

#include "FileHandler.h"

namespace data_manage
{

class MDRFileHandler : public data_manage::FileHandler
{
public:

  MDRFileHandler();

  virtual ~MDRFileHandler();

  /// provides interface for filling dataholder object with genotypes
  void parse_file(std::string filename, Dataholder* holder, int missingValue,
    float statusMissingValue, bool contains_id=false);

  /// parses testing and training files
  void parse_file(std::string train_file, std::string test_file, Dataholder* holder,
	  int missingValue, float statusMissingValue, bool contains_id=false);

  /// provides interface for writing out dataholder object
  void write_file(std::string filename, Dataholder* holder);
  
private:
  int dummy_id;

};

}

#endif /*MDRFILEHANDLER_H_*/
